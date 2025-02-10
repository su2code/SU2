/*!
 * \file CConfig.cpp
 * \brief Main file for managing the config file
 * \author F. Palacios, T. Economon, B. Tracey, H. Kline
 * \version 8.1.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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

#define ENABLE_MAPS
#include <utility>

#include "../include/CConfig.hpp"
#undef ENABLE_MAPS

#include "../include/fem/fem_gauss_jacobi_quadrature.hpp"
#include "../include/fem/fem_geometry_structure.hpp"

#include "../include/basic_types/ad_structure.hpp"
#include "../include/toolboxes/printing_toolbox.hpp"

using namespace PrintingToolbox;

#ifdef PROFILE
#ifdef HAVE_MKL
#include "mkl.h"
#endif
#endif

vector<string> Profile_Function_tp;       /*!< \brief Vector of string names for profiled functions. */
vector<double> Profile_Time_tp;           /*!< \brief Vector of elapsed time for profiled functions. */
vector<double> Profile_ID_tp;             /*!< \brief Vector of group ID number for profiled functions. */
map<string, vector<int> > Profile_Map_tp; /*!< \brief Map containing the final results for profiled functions. */

map<CLong3T, int> GEMM_Profile_MNK;       /*!< \brief Map, which maps the GEMM size to the index where
                                                      the data for this GEMM is stored in several vectors. */
vector<long>   GEMM_Profile_NCalls;       /*!< \brief Vector, which stores the number of calls to this
                                                      GEMM size. */
vector<double> GEMM_Profile_TotTime;      /*!< \brief Total time spent for this GEMM size. */
vector<double> GEMM_Profile_MinTime;      /*!< \brief Minimum time spent for this GEMM size. */
vector<double> GEMM_Profile_MaxTime;      /*!< \brief Maximum time spent for this GEMM size. */

//#pragma omp threadprivate(Profile_Function_tp, Profile_Time_tp, Profile_ID_tp, Profile_Map_tp)


CConfig::CConfig(char case_filename[MAX_STRING_SIZE], SU2_COMPONENT val_software, bool verb_high) {

  /*--- Set the case name to the base config file name without extension ---*/

  caseName = PrintingToolbox::split(string(case_filename),'.')[0];

  base_config = true;

  /*--- Store MPI rank and size ---*/

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  iZone = 0;
  nZone = 1;

  Init();

  /*--- Parsing the config file  ---*/

  SetConfig_Parsing(case_filename);

  /*--- Set the default values for all of the options that weren't set ---*/

  SetDefault();

  /*--- Set number of zone ---*/

  SetnZone();

  /*--- Configuration file postprocessing ---*/

  SetPostprocessing(val_software, iZone, 0);

  /*--- Configuration file boundaries/markers setting ---*/

  SetMarkers(val_software);

  /*--- Configuration file output ---*/

  if ((rank == MASTER_NODE) && verb_high)
    SetOutput(val_software, iZone);

}

CConfig::CConfig(istream &case_buffer, SU2_COMPONENT val_software, bool verb_high) {

  base_config = true;

  iZone = 0;
  nZone = 1;

  Init();

  /*--- Parsing the config file  ---*/

  SetConfig_Parsing(case_buffer);

  /*--- Set the default values for all of the options that weren't set ---*/

  SetDefault();

  /*--- Set number of zone ---*/

  SetnZone();

  /*--- Configuration file postprocessing ---*/

  SetPostprocessing(val_software, iZone, 0);

  /*--- Configuration file boundaries/markers setting ---*/

  SetMarkers(val_software);

  /*--- Configuration file output ---*/

  if ((rank == MASTER_NODE) && verb_high)
    SetOutput(val_software, iZone);

}

CConfig::CConfig(CConfig* config, char case_filename[MAX_STRING_SIZE], SU2_COMPONENT val_software, unsigned short val_iZone, unsigned short val_nZone, bool verb_high) {

  caseName = config->GetCaseName();

  unsigned short val_nDim;

  base_config = false;

  iZone = val_iZone;
  nZone = val_nZone;

  Init();

  /*--- Parsing the config file  ---*/

  SetConfig_Parsing(case_filename);

  /*--- Set default options from base config ---*/

  SetDefaultFromConfig(config);

  /*--- Set the default values for all of the options that weren't set ---*/

  SetDefault();

  /*--- Get the dimension --- */

  val_nDim = GetnDim(Mesh_FileName, Mesh_FileFormat);

  /*--- Configuration file postprocessing ---*/

  SetPostprocessing(val_software, val_iZone, val_nDim);

  /*--- Configuration file boundaries/markers setting ---*/

  SetMarkers(val_software);

  /*--- Configuration file output ---*/

  if ((rank == MASTER_NODE) && verb_high)
    SetOutput(val_software, val_iZone);

  Multizone_Problem = config->GetMultizone_Problem();

}

CConfig::CConfig(char case_filename[MAX_STRING_SIZE], SU2_COMPONENT val_software) {

  /*--- Set the case name to the base config file name without extension ---*/

  caseName = PrintingToolbox::split(string(case_filename),'.')[0];

  base_config = true;

  nZone = 1;
  iZone = 0;

  Init();

  /*--- Parsing the config file  ---*/

  SetConfig_Parsing(case_filename);

  /*--- Set the default values for all of the options that weren't set ---*/

  SetDefault();

  /*--- Set number of zones --- */

  SetnZone();

  /*--- Configuration file postprocessing ---*/

  SetPostprocessing(val_software, 0, 1);

  /*--- Configuration file boundaries/markers setting ---*/

  SetMarkers(val_software);

  /*--- Print the header --- */

  SetHeader(val_software);

}

CConfig::CConfig(char case_filename[MAX_STRING_SIZE], CConfig *config) {

  /*--- Set the case name to the base config file name without extension ---*/

  caseName = PrintingToolbox::split(string(case_filename),'.')[0];

  base_config = true;

  bool runtime_file = false;

  Init();

  /*--- Parsing the config file  ---*/

  runtime_file = SetRunTime_Parsing(case_filename);

  /*--- Set the default values for all of the options that weren't set ---*/

  SetDefault();

  /*--- Update original config file ---*/

  if (runtime_file) {
    if (all_options.find("TIME_ITER") == all_options.end())
      config->SetnTime_Iter(nTimeIter);
  }
}

SU2_MPI::Comm CConfig::GetMPICommunicator() const {

  return SU2_Communicator;

}

void CConfig::Init(){

  /*--- Store MPI rank and size ---*/

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  /*--- Initialize pointers to Null---*/

  SetPointersNull();

  /*--- Reading config options  ---*/

  SetConfig_Options();

}

void CConfig::SetMPICommunicator(SU2_MPI::Comm Communicator) {

  SU2_Communicator = Communicator;

}

void CConfig::addDoubleOption(const string& name, su2double & option_field, su2double default_value) {
  // Check if the key is already in the map. If this fails, it is coder error
  // and not user error, so throw.
  assert(option_map.find(name) == option_map.end());

  // Add this option to the list of all the options
  all_options.insert(pair<string, bool>(name, true));

  // Create the parser for a su2double option with a reference to the option_field and the desired
  // default value. This will take the string in the config file, convert it to a su2double, and
  // place that su2double in the memory location specified by the reference.
  COptionBase* val = new COptionDouble(name, option_field, default_value);

  // Create an association between the option name ("CFL") and the parser generated above.
  // During configuration, the parsing script will get the option name, and use this map
  // to find how to parse that option.
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addStringOption(const string& name, string & option_field, string default_value) {

  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionString(name, option_field, std::move(default_value));
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addIntegerOption(const string& name, int & option_field, int default_value) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionInt(name, option_field, default_value);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addUnsignedLongOption(const string& name, unsigned long & option_field, unsigned long default_value) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionULong(name, option_field, default_value);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addUnsignedShortOption(const string& name, unsigned short & option_field, unsigned short default_value) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionUShort(name, option_field, default_value);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addLongOption(const string& name, long & option_field, long default_value) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionLong(name, option_field, default_value);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addBoolOption(const string& name, bool & option_field, bool default_value) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionBool(name, option_field, default_value);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

// enum types work differently than all of the others because there are a small number of valid
// string entries for the type. One must also provide a list of all the valid strings of that type.
template <class Tenum, class TField>
void CConfig::addEnumOption(const string name, TField& option_field, const map<string,Tenum>& enum_map, Tenum default_value) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionEnum<Tenum, TField>(name, enum_map, option_field, default_value);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

// input_size is the number of options read in from the config file
template <class Tenum, class TField>
void CConfig::addEnumListOption(const string name, unsigned short& input_size, TField*& option_field, const map<string, Tenum>& enum_map) {
  input_size = 0;
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionEnumList<Tenum,TField>(name, enum_map, option_field, input_size);
  option_map.insert( pair<string, COptionBase*>(name, val) );
}

void CConfig::addDoubleArrayOption(const string& name, const int size, su2double* option_field) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionArray<su2double>(name, size, option_field);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addUShortArrayOption(const string& name, const int size, unsigned short* option_field) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionArray<unsigned short>(name, size, option_field);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addDoubleListOption(const string& name, unsigned short & size, su2double * & option_field) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionDoubleList(name, size, option_field);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addShortListOption(const string& name, unsigned short & size, short * & option_field) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionShortList(name, size, option_field);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addUShortListOption(const string& name, unsigned short & size, unsigned short * & option_field) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionUShortList(name, size, option_field);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addULongListOption(const string& name, unsigned short & size, unsigned long * & option_field) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionULongList(name, size, option_field);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addStringListOption(const string& name, unsigned short & num_marker, string* & option_field) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionStringList(name, num_marker, option_field);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addConvectOption(const string& name, unsigned short & space_field, CENTERED & centered_field, UPWIND & upwind_field) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionConvect(name, space_field, centered_field, upwind_field);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addConvectFEMOption(const string& name, unsigned short & space_field, unsigned short & fem_field) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionFEMConvect(name, space_field, fem_field);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addMathProblemOption(const string& name, bool & ContinuousAdjoint, const bool & ContinuousAdjoint_default,
                          bool & DiscreteAdjoint, const bool & DiscreteAdjoint_default,
                          bool & Restart_Flow, const bool & Restart_Flow_default) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionMathProblem(name, ContinuousAdjoint, ContinuousAdjoint_default, DiscreteAdjoint, DiscreteAdjoint_default, Restart_Flow, Restart_Flow_default);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addDVParamOption(const string& name, unsigned short & nDV_field, su2double** & paramDV, string* & FFDTag,
                      unsigned short* & design_variable) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionDVParam(name, nDV_field, paramDV, FFDTag, design_variable);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addDVValueOption(const string& name, unsigned short* & nDVValue_field, su2double** & valueDV, unsigned short & nDV_field,  su2double** & paramDV,
                      unsigned short* & design_variable) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionDVValue(name, nDVValue_field, valueDV, nDV_field, paramDV, design_variable);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addFFDDefOption(const string& name, unsigned short & nFFD_field, su2double** & coordFFD, string* & FFDTag) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionFFDDef(name, nFFD_field, coordFFD, FFDTag);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addFFDDegreeOption(const string& name, unsigned short & nFFD_field, unsigned short** & degreeFFD) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionFFDDegree(name, nFFD_field, degreeFFD);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addStringDoubleListOption(const string& name, unsigned short & list_size, string * & string_field,
                                        su2double* & double_field) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionStringValuesList<su2double>(name, list_size, string_field, double_field);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addInletOption(const string& name, unsigned short & nMarker_Inlet, string * & Marker_Inlet,
                    su2double* & Ttotal, su2double* & Ptotal, su2double** & FlowDir) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionInlet(name, nMarker_Inlet, Marker_Inlet, Ttotal, Ptotal, FlowDir);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addInletSpeciesOption(const string& name, unsigned short & nMarker_Inlet_Species,
                                    string * & Marker_Inlet_Species, su2double** & inlet_species_val,
                                    unsigned short & nSpecies_per_Inlet) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionStringValuesList<su2double*>(name, nMarker_Inlet_Species, Marker_Inlet_Species,
                                                             inlet_species_val, nSpecies_per_Inlet);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addInletTurbOption(const string& name, unsigned short& nMarker_Inlet_Turb, string*& Marker_Inlet_Turb,
                                 su2double**& Turb_Properties_val, unsigned short& nTurb_Properties) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionStringValuesList<su2double*>(name, nMarker_Inlet_Turb, Marker_Inlet_Turb,
                                                             Turb_Properties_val, nTurb_Properties);
  option_map.insert(pair<string, COptionBase*>(name, val));
}

template <class Tenum>
void CConfig::addRiemannOption(const string name, unsigned short & nMarker_Riemann, string * & Marker_Riemann, unsigned short* & option_field, const map<string, Tenum> & enum_map,
                               su2double* & var1, su2double* & var2, su2double** & FlowDir) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionRiemann<Tenum>(name, nMarker_Riemann, Marker_Riemann, option_field, enum_map, var1, var2, FlowDir);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

template <class Tenum>
void CConfig::addGilesOption(const string name, unsigned short & nMarker_Giles, string * & Marker_Giles, unsigned short* & option_field, const map<string, Tenum> & enum_map,
                             su2double* & var1, su2double* & var2, su2double** & FlowDir, su2double* & relaxfactor1, su2double* & relaxfactor2) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionGiles<Tenum>(name, nMarker_Giles, Marker_Giles, option_field, enum_map, var1, var2, FlowDir, relaxfactor1, relaxfactor2);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addExhaustOption(const string& name, unsigned short & nMarker_Exhaust, string * & Marker_Exhaust,
                               su2double* & Ttotal, su2double* & Ptotal) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionExhaust(name, nMarker_Exhaust, Marker_Exhaust, Ttotal, Ptotal);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addPeriodicOption(const string & name, unsigned short & nMarker_PerBound,
                                string* & Marker_PerBound, string* & Marker_PerDonor,
                                su2double** & RotCenter, su2double** & RotAngles, su2double** & Translation) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionPeriodic(name, nMarker_PerBound, Marker_PerBound, Marker_PerDonor, RotCenter, RotAngles, Translation);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addTurboPerfOption(const string & name, unsigned short & nMarker_TurboPerf,
                                 string* & Marker_TurboBoundIn, string* & Marker_TurboBoundOut, string* & Marker_Turbomachinery) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionTurboPerformance(name, nMarker_TurboPerf, Marker_TurboBoundIn, Marker_TurboBoundOut, Marker_Turbomachinery);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addActDiskOption(const string & name, unsigned short & nMarker_ActDiskInlet,
                               unsigned short & nMarker_ActDiskOutlet, string* & Marker_ActDiskInlet,
                               string* & Marker_ActDiskOutlet, su2double** & ActDisk_PressJump,
                               su2double** & ActDisk_TempJump, su2double** & ActDisk_Omega) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionActDisk(name, nMarker_ActDiskInlet, nMarker_ActDiskOutlet, Marker_ActDiskInlet,
                                        Marker_ActDiskOutlet, ActDisk_PressJump, ActDisk_TempJump, ActDisk_Omega);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addActDiskBemOption(const string& name,
                                  unsigned short& nMarker_ActDiskBemInlet, unsigned short& nMarker_ActDiskBemOutlet,
                                  string*& Marker_ActDiskBemInlet, string*& Marker_ActDiskBemOutlet,
                                  su2double**& ActDiskBem_X, su2double**& ActDiskBem_Y, su2double**& ActDiskBem_Z) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionActDisk(name,
                                        nMarker_ActDiskBemInlet, nMarker_ActDiskBemOutlet,
                                        Marker_ActDiskBemInlet, Marker_ActDiskBemOutlet,
                                        ActDiskBem_X, ActDiskBem_Y, ActDiskBem_Z);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addWallFunctionOption(const string &name, unsigned short &list_size, string* &string_field,
                                    WALL_FUNCTIONS* &val_Kind_WF, unsigned short** &val_IntInfo_WF,
                                    su2double** &val_DoubleInfo_WF) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionWallFunction(name, list_size, string_field, val_Kind_WF,
                                             val_IntInfo_WF, val_DoubleInfo_WF);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

void CConfig::addPythonOption(const string& name) {
  assert(option_map.find(name) == option_map.end());
  all_options.insert(pair<string, bool>(name, true));
  COptionBase* val = new COptionPython(name);
  option_map.insert(pair<string, COptionBase *>(name, val));
}

unsigned short CConfig::GetnZone(const string& val_mesh_filename, unsigned short val_format) {

  int nZone = 1; /* Default value if nothing is specified. */

  switch (val_format) {
    case SU2: {

      /*--- Local variables for reading the SU2 file. ---*/
      string text_line;
      ifstream mesh_file;

      /*--- Check if the mesh file can be opened for reading. ---*/
      mesh_file.open(val_mesh_filename.c_str(), ios::in);
      if (mesh_file.fail())
        SU2_MPI::Error(string("There is no geometry file called ") + val_mesh_filename,
                              CURRENT_FUNCTION);

      /*--- Read the SU2 mesh file until the zone data is reached or
            when it can be decided that it is not present. ---*/
      while( getline (mesh_file, text_line) ) {

        /*--- Search for the "NZONE" keyword to see if there are multiple Zones ---*/
        if(text_line.find ("NZONE=",0) != string::npos) {
          text_line.erase (0,6); nZone = atoi(text_line.c_str());
          break;
        }

        /*--- If one of the keywords IZONE, NELEM or NPOIN, NMARK is encountered,
              it can be assumed that the NZONE keyword is not present and the loop
              can be terminated. ---*/
        if(text_line.find ("IZONE=",0) != string::npos) break;
        if(text_line.find ("NELEM=",0) != string::npos) break;
        if(text_line.find ("NPOIN=",0) != string::npos) break;
        if(text_line.find ("NMARK=",0) != string::npos) break;
      }

      mesh_file.close();
      break;

    }

    case CGNS_GRID: {

#ifdef HAVE_CGNS

      /*--- Local variables which are needed when calling the CGNS mid-level API. ---*/

      int fn, nbases = 0, nzones = 0, file_type;
      int cell_dim = 0, phys_dim = 0;
      char basename[CGNS_STRING_SIZE];

      /*--- Check whether the supplied file is truly a CGNS file. ---*/

      if ( cg_is_cgns(val_mesh_filename.c_str(), &file_type) != CG_OK ) {
        SU2_MPI::Error(val_mesh_filename +
                       string(" was not found or is not a properly formatted CGNS file.\n") +
                       string("Note that SU2 expects unstructured CGNS files in ADF data format."),
                       CURRENT_FUNCTION);
      }

      /*--- Open the CGNS file for reading. The value of fn returned
       is the specific index number for this file and will be
       repeatedly used in the function calls. ---*/

      if (cg_open(val_mesh_filename.c_str(), CG_MODE_READ, &fn)) cg_error_exit();

      /*--- Get the number of databases. This is the highest node
       in the CGNS heirarchy. ---*/

      if (cg_nbases(fn, &nbases)) cg_error_exit();

      /*--- Check if there is more than one database. Throw an
       error if there is because this reader can currently
       only handle one database. ---*/

      if ( nbases > 1 ) {
        SU2_MPI::Error("CGNS reader currently incapable of handling more than 1 database." ,
                       CURRENT_FUNCTION);
      }

      /*--- Read the databases. Note that the indexing starts at 1. ---*/

      for ( int i = 1; i <= nbases; i++ ) {

        if (cg_base_read(fn, i, basename, &cell_dim, &phys_dim)) cg_error_exit();

        /*--- Get the number of zones for this base. ---*/

        if (cg_nzones(fn, i, &nzones)) cg_error_exit();

      }

      /*--- Close the CGNS file. ---*/

      if ( cg_close(fn) ) cg_error_exit();

      /*--- Set the number of zones as read from the CGNS file ---*/

      nZone = nzones;

#else
      SU2_MPI::Error(string(" SU2 built without CGNS support. \n") +
                     string(" To use CGNS, build SU2 accordingly."),
                     CURRENT_FUNCTION);
#endif

      break;
    }
    case RECTANGLE: {
      nZone = 1;
      break;
    }
    case BOX: {
      nZone = 1;
      break;
    }
  }

  return (unsigned short) nZone;

}

unsigned short CConfig::GetnDim(const string& val_mesh_filename, unsigned short val_format) {

  short nDim = -1;

  switch (val_format) {
    case SU2: {

      /*--- Local variables for reading the SU2 file. ---*/
      string text_line;
      ifstream mesh_file;

      /*--- Open grid file ---*/
      mesh_file.open(val_mesh_filename.c_str(), ios::in);
      if (mesh_file.fail()) {
        SU2_MPI::Error(string("The SU2 mesh file named ") + val_mesh_filename + string(" was not found."), CURRENT_FUNCTION);
      }

      /*--- Read the SU2 mesh file until the dimension data is reached
            or when it can be decided that it is not present. ---*/
      while( getline (mesh_file, text_line) ) {

        /*--- Search for the "NDIME" keyword to determine the number
              of dimensions.  ---*/
        if(text_line.find ("NDIME=",0) != string::npos) {
          text_line.erase (0,6); nDim = atoi(text_line.c_str());
          break;
        }

        /*--- If one of the keywords NELEM or NPOIN, NMARK is encountered,
              it can be assumed that the NZONE keyword is not present and
              the loop can be terminated. ---*/
        if(text_line.find ("NELEM=",0) != string::npos) break;
        if(text_line.find ("NPOIN=",0) != string::npos) break;
        if(text_line.find ("NMARK=",0) != string::npos) break;
      }

      mesh_file.close();

      /*--- Throw an error if the dimension was not found. ---*/
      if (nDim == -1) {
        SU2_MPI::Error(val_mesh_filename + string(" is not an SU2 mesh file or has the wrong format \n ('NDIME=' not found). Please check."),
                       CURRENT_FUNCTION);
      }

      break;
    }

    case CGNS_GRID: {

#ifdef HAVE_CGNS

      /*--- Local variables which are needed when calling the CGNS mid-level API. ---*/
      int fn, nbases, file_type;
      int cell_dim, phys_dim;
      char basename[CGNS_STRING_SIZE];

      /*--- Check whether the supplied file is truly a CGNS file. ---*/
      if ( cg_is_cgns(val_mesh_filename.c_str(), &file_type) != CG_OK ) {
        SU2_MPI::Error(val_mesh_filename +
                       string(" was not found or is not a properly formatted CGNS file.\n") +
                       string("Note that SU2 expects unstructured CGNS files in ADF data format."),
                       CURRENT_FUNCTION);
      }

      /*--- Open the CGNS file for reading. The value of fn returned
            is the specific index number for this file and will be
            repeatedly used in the function calls. ---*/
      if (cg_open(val_mesh_filename.c_str(), CG_MODE_READ, &fn) != CG_OK) cg_error_exit();

      /*--- Get the number of databases. This is the highest node
            in the CGNS heirarchy. ---*/
      if (cg_nbases(fn, &nbases) != CG_OK) cg_error_exit();

      /*--- Check if there is more than one database. Throw an
            error if there is because this reader can currently
            only handle one database. ---*/
      if ( nbases > 1 )
        SU2_MPI::Error("CGNS reader currently incapable of handling more than 1 database." ,
                       CURRENT_FUNCTION);

      /*--- Read the database. Note that the indexing starts at 1.
            Afterwards close the file again. ---*/
      if (cg_base_read(fn, 1, basename, &cell_dim, &phys_dim) != CG_OK) cg_error_exit();
      if (cg_close(fn) != CG_OK) cg_error_exit();

      /*--- Set the problem dimension as read from the CGNS file ---*/
      nDim = cell_dim;

#else
      SU2_MPI::Error(string(" SU2 built without CGNS support. \n") +
                     string(" To use CGNS, build SU2 accordingly."),
                     CURRENT_FUNCTION);
#endif

      break;
    }
    case RECTANGLE: {
      nDim = 2;
      break;
    }
    case BOX: {
      nDim = 3;
      break;
    }
  }

  /*--- After reading the mesh, assert that the dimension is equal to 2 or 3. ---*/
  assert((nDim == 2) || (nDim == 3));

  return (unsigned short) nDim;
}

void CConfig::SetPointersNull() {

  Marker_CfgFile_GeoEval      = nullptr;   Marker_All_GeoEval       = nullptr;
  Marker_CfgFile_Monitoring   = nullptr;   Marker_All_Monitoring    = nullptr;
  Marker_CfgFile_Designing    = nullptr;   Marker_All_Designing     = nullptr;
  Marker_CfgFile_Plotting     = nullptr;   Marker_All_Plotting      = nullptr;
  Marker_CfgFile_Analyze      = nullptr;   Marker_All_Analyze       = nullptr;
  Marker_CfgFile_DV           = nullptr;   Marker_All_DV            = nullptr;
  Marker_CfgFile_Moving       = nullptr;   Marker_All_Moving        = nullptr;
  Marker_CfgFile_PerBound     = nullptr;   Marker_All_PerBound      = nullptr;    Marker_PerBound   = nullptr;
  Marker_CfgFile_Turbomachinery = nullptr; Marker_All_Turbomachinery = nullptr;
  Marker_CfgFile_TurbomachineryFlag = nullptr; Marker_All_TurbomachineryFlag = nullptr;
  Marker_CfgFile_MixingPlaneInterface = nullptr; Marker_All_MixingPlaneInterface = nullptr;
  Marker_CfgFile_ZoneInterface = nullptr;
  Marker_CfgFile_Deform_Mesh   = nullptr;  Marker_All_Deform_Mesh   = nullptr;
  Marker_CfgFile_Deform_Mesh_Sym_Plane   = nullptr;  Marker_All_Deform_Mesh_Sym_Plane   = nullptr;
  Marker_CfgFile_Fluid_Load    = nullptr;  Marker_All_Fluid_Load    = nullptr;
  Marker_CfgFile_SobolevBC     = nullptr;  Marker_All_SobolevBC     = nullptr;

  Marker_CfgFile_Turbomachinery       = nullptr; Marker_All_Turbomachinery       = nullptr;
  Marker_CfgFile_TurbomachineryFlag   = nullptr; Marker_All_TurbomachineryFlag   = nullptr;
  Marker_CfgFile_MixingPlaneInterface = nullptr; Marker_All_MixingPlaneInterface = nullptr;

  Marker_CfgFile_PyCustom     = nullptr;   Marker_All_PyCustom      = nullptr;

  Marker_DV                   = nullptr;   Marker_Moving            = nullptr;    Marker_Monitoring = nullptr;
  Marker_Designing            = nullptr;   Marker_GeoEval           = nullptr;    Marker_Plotting   = nullptr;
  Marker_Analyze              = nullptr;   Marker_PyCustom          = nullptr;    Marker_WallFunctions        = nullptr;
  Marker_CfgFile_KindBC       = nullptr;   Marker_All_KindBC        = nullptr;    Marker_SobolevBC  = nullptr;
  Marker_StrongBC             = nullptr;

  Kind_WallFunctions       = nullptr;
  IntInfo_WallFunctions    = nullptr;
  DoubleInfo_WallFunctions = nullptr;

  Config_Filenames = nullptr;

  /*--- Marker Pointers ---*/

  Marker_Euler                = nullptr;    Marker_FarField             = nullptr;    Marker_Custom              = nullptr;
  Marker_SymWall              = nullptr;    Marker_PerBound             = nullptr;
  Marker_PerDonor             = nullptr;    Marker_NearFieldBound       = nullptr;    Marker_Inlet_Turb          = nullptr;
  Marker_Deform_Mesh          = nullptr;    Marker_Deform_Mesh_Sym_Plane= nullptr;    Marker_Fluid_Load          = nullptr;
  Marker_Inlet                = nullptr;    Marker_Outlet               = nullptr;    Marker_Inlet_Species       = nullptr;
  Marker_Supersonic_Inlet     = nullptr;    Marker_Supersonic_Outlet    = nullptr;    Marker_Smoluchowski_Maxwell= nullptr;
  Marker_Isothermal           = nullptr;    Marker_HeatFlux             = nullptr;    Marker_EngineInflow        = nullptr;
  Marker_Load                 = nullptr;    Marker_Disp_Dir             = nullptr;    Marker_RoughWall           = nullptr;
  Marker_EngineExhaust        = nullptr;    Marker_Displacement         = nullptr;    Marker_Load                = nullptr;
  Marker_Load_Dir             = nullptr;    Marker_Clamped             = nullptr;
  Marker_Internal             = nullptr;
  Marker_All_TagBound         = nullptr;    Marker_CfgFile_TagBound     = nullptr;    Marker_All_KindBC          = nullptr;
  Marker_CfgFile_KindBC       = nullptr;    Marker_All_SendRecv         = nullptr;    Marker_All_PerBound        = nullptr;
  Marker_ZoneInterface        = nullptr;    Marker_All_ZoneInterface    = nullptr;    Marker_Riemann             = nullptr;
  Marker_Fluid_InterfaceBound = nullptr;    Marker_CHTInterface         = nullptr;    Marker_Damper              = nullptr;
  Marker_Emissivity           = nullptr;    Marker_HeatTransfer         = nullptr;

    /*--- Boundary Condition settings ---*/

  Isothermal_Temperature = nullptr;    HeatTransfer_Coeff     = nullptr;    HeatTransfer_WallTemp  = nullptr;
  Heat_Flux              = nullptr;    Displ_Value            = nullptr;    Load_Value             = nullptr;
  Damper_Constant        = nullptr;    Wall_Emissivity        = nullptr;
  Roughness_Height       = nullptr;

  /*--- Inlet Outlet Boundary Condition settings ---*/

  Inlet_Ttotal    = nullptr;    Inlet_Ptotal      = nullptr;
  Inlet_FlowDir   = nullptr;    Inlet_Temperature = nullptr;    Inlet_Pressure = nullptr;
  Inlet_Velocity  = nullptr;
  Outlet_Pressure = nullptr;    Inlet_SpeciesVal  = nullptr;    Inlet_TurbVal = nullptr;

  /*--- Engine Boundary Condition settings ---*/

  Inflow_Pressure      = nullptr;    Inflow_MassFlow    = nullptr;    Inflow_ReverseMassFlow  = nullptr;
  Inflow_TotalPressure = nullptr;    Inflow_Temperature = nullptr;    Inflow_TotalTemperature = nullptr;
  Inflow_RamDrag       = nullptr;    Inflow_Force       = nullptr;    Inflow_Power            = nullptr;
  Inflow_Mach          = nullptr;

  Exhaust_Pressure        = nullptr;   Exhaust_Temperature        = nullptr;    Exhaust_MassFlow = nullptr;
  Exhaust_TotalPressure   = nullptr;   Exhaust_TotalTemperature   = nullptr;
  Exhaust_GrossThrust     = nullptr;   Exhaust_Force              = nullptr;
  Exhaust_Power           = nullptr;   Exhaust_Temperature_Target = nullptr;
  Exhaust_Pressure_Target = nullptr;

  Engine_Mach  = nullptr;    Engine_Force        = nullptr;
  Engine_Power = nullptr;    Engine_NetThrust    = nullptr;    Engine_GrossThrust = nullptr;
  Engine_Area  = nullptr;    EngineInflow_Target = nullptr;

  Exhaust_Temperature_Target  = nullptr;   Exhaust_Temperature     = nullptr;   Exhaust_Pressure      = nullptr;
  Exhaust_Pressure_Target     = nullptr;   Inlet_Ttotal            = nullptr;   Inlet_Ptotal          = nullptr;
  Inlet_FlowDir               = nullptr;   Inlet_Temperature       = nullptr;   Inlet_Pressure        = nullptr;
  Inlet_Velocity              = nullptr;   Inflow_Mach             = nullptr;   Inflow_Pressure       = nullptr;
  Outlet_Pressure             = nullptr;   Isothermal_Temperature  = nullptr;

  ElasticityMod = nullptr;
  PoissonRatio = nullptr;
  MaterialDensity = nullptr;
  MaterialThermalExpansion = nullptr;

  Load_Dir = nullptr;            Load_Dir_Value = nullptr;          Load_Dir_Multiplier = nullptr;
  Disp_Dir = nullptr;            Disp_Dir_Value = nullptr;          Disp_Dir_Multiplier = nullptr;
  Electric_Field_Mod = nullptr;  Electric_Field_Dir = nullptr;      RefNode_Displacement = nullptr;

  Electric_Constant = nullptr;

  /*--- Actuator Disk Boundary Condition settings ---*/

  ActDiskInlet_Pressure         = nullptr;    ActDiskInlet_TotalPressure = nullptr;    ActDiskInlet_Temperature = nullptr;
  ActDiskInlet_TotalTemperature = nullptr;    ActDiskInlet_MassFlow      = nullptr;    ActDiskInlet_RamDrag     = nullptr;
  ActDiskInlet_Force            = nullptr;    ActDiskInlet_Power         = nullptr;

  ActDiskOutlet_Pressure      = nullptr;
  ActDiskOutlet_TotalPressure = nullptr;   ActDiskOutlet_GrossThrust = nullptr;  ActDiskOutlet_Force            = nullptr;
  ActDiskOutlet_Power         = nullptr;   ActDiskOutlet_Temperature = nullptr;  ActDiskOutlet_TotalTemperature = nullptr;
  ActDiskOutlet_MassFlow      = nullptr;

  ActDiskOutlet_Thrust_BEM = nullptr;
  ActDiskOutlet_Torque_BEM = nullptr;

  ActDisk_DeltaPress      = nullptr;    ActDisk_DeltaTemp      = nullptr;
  ActDisk_TotalPressRatio = nullptr;    ActDisk_TotalTempRatio = nullptr;    ActDisk_StaticPressRatio = nullptr;
  ActDisk_StaticTempRatio = nullptr;    ActDisk_NetThrust      = nullptr;    ActDisk_GrossThrust      = nullptr;
  ActDisk_Power           = nullptr;    ActDisk_MassFlow       = nullptr;    ActDisk_Area             = nullptr;
  ActDisk_ReverseMassFlow = nullptr;    Surface_MassFlow        = nullptr;   Surface_Mach             = nullptr;
  Surface_Temperature      = nullptr;   Surface_Pressure         = nullptr;  Surface_Density          = nullptr;   Surface_Enthalpy          = nullptr;
  Surface_NormalVelocity   = nullptr;   Surface_TotalTemperature = nullptr;  Surface_TotalPressure    = nullptr;   Surface_PressureDrop    = nullptr;
  Surface_DC60             = nullptr;   Surface_IDC = nullptr;
  Surface_Species_Variance = nullptr;   Surface_Species_0 = nullptr;

  Outlet_MassFlow      = nullptr;       Outlet_Density      = nullptr;      Outlet_Area     = nullptr;

  Surface_Uniformity = nullptr; Surface_SecondaryStrength = nullptr; Surface_SecondOverUniform = nullptr;
  Surface_MomentumDistortion = nullptr;

  Surface_IDC_Mach        = nullptr;    Surface_IDR            = nullptr;    ActDisk_Mach             = nullptr;
  ActDisk_Force           = nullptr;    ActDisk_BCThrust       = nullptr;    ActDisk_BCThrust_Old     = nullptr;

  /*--- Miscellaneous/unsorted ---*/

  Aeroelastic_plunge  = nullptr;
  Aeroelastic_pitch   = nullptr;

  CFL_AdaptParam      = nullptr;
  CFL                 = nullptr;
  PlaneTag            = nullptr;
  ParamDV             = nullptr;
  DV_Value            = nullptr;
  Design_Variable     = nullptr;

  TimeDOFsADER_DG           = nullptr;
  TimeIntegrationADER_DG    = nullptr;
  WeightsIntegrationADER_DG = nullptr;
  RK_Alpha_Step             = nullptr;
  MG_CorrecSmooth           = nullptr;
  MG_PreSmooth              = nullptr;
  MG_PostSmooth             = nullptr;
  Int_Coeffs                = nullptr;

  Kind_Inc_Inlet = nullptr;
  Kind_Inc_Outlet = nullptr;

  Kind_ObjFunc   = nullptr;

  Weight_ObjFunc = nullptr;

  /*--- Species solver pointers. ---*/

  Species_Init           = nullptr;
  Species_Clipping_Min   = nullptr;
  Species_Clipping_Max   = nullptr;
  spark_reaction_rates   = nullptr;

  /*--- Moving mesh pointers ---*/

  nKind_SurfaceMovement = 0;
  Kind_SurfaceMovement = nullptr;
  LocationStations   = nullptr;
  MarkerMotion_Origin     = nullptr;
  MarkerTranslation_Rate  = nullptr;
  MarkerRotation_Rate     = nullptr;
  MarkerPitching_Omega    = nullptr;
  MarkerPitching_Ampl     = nullptr;
  MarkerPitching_Phase    = nullptr;
  MarkerPlunging_Omega    = nullptr;
  MarkerPlunging_Ampl     = nullptr;
  RefOriginMoment_X   = nullptr;    RefOriginMoment_Y   = nullptr;    RefOriginMoment_Z   = nullptr;
  MoveMotion_Origin   = nullptr;

  /*--- Periodic BC pointers. ---*/

  Periodic_Translation= nullptr;    Periodic_RotAngles  = nullptr;    Periodic_RotCenter  = nullptr;

  /* Harmonic Balance Frequency pointer */

  Omega_HB = nullptr;

  /*--- Initialize some default arrays to NULL. ---*/

  Riemann_FlowDir       = nullptr;
  Giles_FlowDir         = nullptr;
  CoordFFDBox           = nullptr;
  DegreeFFDBox          = nullptr;
  FFDTag                = nullptr;
  nDV_Value             = nullptr;
  TagFFDBox             = nullptr;

  Kind_Data_Riemann        = nullptr;
  Riemann_Var1             = nullptr;
  Riemann_Var2             = nullptr;
  Kind_Data_Giles          = nullptr;
  Giles_Var1               = nullptr;
  Giles_Var2               = nullptr;
  RelaxFactorAverage       = nullptr;
  RelaxFactorFourier       = nullptr;
  nSpan_iZones             = nullptr;
  Kind_TurboMachinery      = nullptr;

  Marker_MixingPlaneInterface  = nullptr;
  Marker_TurboBoundIn          = nullptr;
  Marker_TurboBoundOut         = nullptr;
  Marker_Turbomachinery        = nullptr;
  Marker_Giles                 = nullptr;
  Marker_Shroud                = nullptr;

  nBlades                      = nullptr;
  FreeStreamTurboNormal        = nullptr;

  top_optim_kernels       = nullptr;
  top_optim_kernel_params = nullptr;
  top_optim_filter_radius = nullptr;

  ScreenOutput = nullptr;
  HistoryOutput = nullptr;
  VolumeOutput = nullptr;
  VolumeOutputFiles = nullptr;
  VolumeOutputFrequencies = nullptr;
  ConvField = nullptr;

  /*--- Variable initialization ---*/

  TimeIter   = 0;
  InnerIter  = 0;
  nIntCoeffs = 0;
  OuterIter  = 0;

  AoA_Offset = 0;
  AoS_Offset = 0;

  nMarker_PerBound = 0;

  Aeroelastic_Simulation = false;

  nSpanMaxAllZones = 1;

  Restart_Bandwidth_Agg = 0.0;

  Mesh_Box_Size = nullptr;

  Time_Ref = 1.0;

  Delta_UnstTime = 0.0;
  Delta_UnstTimeND = 0.0;
  Total_UnstTime = 0.0;
  Total_UnstTimeND = 0.0;

  Kind_TimeNumScheme = EULER_IMPLICIT;

}

void CConfig::SetConfig_Options() {

  // This config file is parsed by a number of programs to make it easy to write SU2
  // wrapper scripts (in python, go, etc.) so please do
  // the best you can to follow the established format. It's very hard to parse c++ code
  // and none of us that write the parsers want to write a full c++ interpreter. Please
  // play nice with the existing format so that you don't break the existing scripts.

  /* BEGIN_CONFIG_OPTIONS */

  /*!\par CONFIG_CATEGORY: Problem Definition \ingroup Config */
  /*--- Options related to problem definition and partitioning ---*/

  /*!\brief SOLVER \n DESCRIPTION: Type of solver \n Options: see \link Solver_Map \endlink \n DEFAULT: NONE \ingroup Config*/
  addEnumOption("SOLVER", Kind_Solver, Solver_Map, MAIN_SOLVER::NONE);
  /*!\brief MULTIZONE \n DESCRIPTION: Enable multizone mode \ingroup Config*/
  addBoolOption("MULTIZONE", Multizone_Problem, NO);
  /*!\brief PHYSICAL_PROBLEM \n DESCRIPTION: Physical governing equations \n Options: see \link Solver_Map \endlink \n DEFAULT: NONE \ingroup Config*/
  addEnumOption("MULTIZONE_SOLVER", Kind_MZSolver, Multizone_Map, ENUM_MULTIZONE::MZ_BLOCK_GAUSS_SEIDEL);
#ifdef CODI_REVERSE_TYPE
  const bool discAdjDefault = true;
#else
  const bool discAdjDefault = false;
#endif
  /*!\brief MATH_PROBLEM  \n DESCRIPTION: Mathematical problem \n  Options: DIRECT, ADJOINT \ingroup Config*/
  addMathProblemOption("MATH_PROBLEM", ContinuousAdjoint, false, DiscreteAdjoint, discAdjDefault, Restart_Flow, discAdjDefault);
  /*!\brief KIND_TURB_MODEL \n DESCRIPTION: Specify turbulence model \n Options: see \link Turb_Model_Map \endlink \n DEFAULT: NONE \ingroup Config*/
  addEnumOption("KIND_TURB_MODEL", Kind_Turb_Model, Turb_Model_Map, TURB_MODEL::NONE);
  /*!\brief SST_OPTIONS \n DESCRIPTION: Specify SST turbulence model options/corrections. \n Options: see \link SST_Options_Map \endlink \n DEFAULT: NONE \ingroup Config*/
  addEnumListOption("SST_OPTIONS", nSST_Options, SST_Options, SST_Options_Map);
  /*!\brief SST_OPTIONS \n DESCRIPTION: Specify SA turbulence model options/corrections. \n Options: see \link SA_Options_Map \endlink \n DEFAULT: NONE \ingroup Config*/
  addEnumListOption("SA_OPTIONS", nSA_Options, SA_Options, SA_Options_Map);

  /*!\brief KIND_TRANS_MODEL \n DESCRIPTION: Specify transition model OPTIONS: see \link Trans_Model_Map \endlink \n DEFAULT: NONE \ingroup Config*/
  addEnumOption("KIND_TRANS_MODEL", Kind_Trans_Model, Trans_Model_Map, TURB_TRANS_MODEL::NONE);
  /*!\brief SST_OPTIONS \n DESCRIPTION: Specify LM transition model options/correlations. \n Options: see \link LM_Options_Map \endlink \n DEFAULT: NONE \ingroup Config*/
  addEnumListOption("LM_OPTIONS", nLM_Options, LM_Options, LM_Options_Map);
  /*!\brief AFT_OPTIONS \n DESCRIPTION: Specify AFT transition model options/correlations. \n Options: see \link AFT_Options_Map \endlink \n DEFAULT: NONE \ingroup Config*/
  addEnumListOption("AFT_OPTIONS", nAFT_Options, AFT_Options, AFT_Options_Map);
  /*!\brief HROUGHNESS \n DESCRIPTION: Value of RMS roughness for transition model \n DEFAULT: 1E-6 \ingroup Config*/
  addDoubleOption("HROUGHNESS", hRoughness, 1e-6);

  /*!\brief KIND_SCALAR_MODEL \n DESCRIPTION: Specify scalar transport model \n Options: see \link Scalar_Model_Map \endlink \n DEFAULT: NONE \ingroup Config*/
  addEnumOption("KIND_SCALAR_MODEL", Kind_Species_Model, Species_Model_Map, SPECIES_MODEL::NONE);

  /*!\brief KIND_SGS_MODEL \n DESCRIPTION: Specify subgrid scale model OPTIONS: see \link SGS_Model_Map \endlink \n DEFAULT: NONE \ingroup Config*/
  addEnumOption("KIND_SGS_MODEL", Kind_SGS_Model, SGS_Model_Map, TURB_SGS_MODEL::NONE);

  /*!\brief KIND_FEM_DG_SHOCK \n DESCRIPTION: Specify shock capturing method for DG OPTIONS: see \link ShockCapturingDG_Map \endlink \n DEFAULT: NONE \ingroup Config*/
  addEnumOption("KIND_FEM_DG_SHOCK", Kind_FEM_Shock_Capturing_DG, ShockCapturingDG_Map, FEM_SHOCK_CAPTURING_DG::NONE);

  /*!\brief KIND_VERIFICATION_SOLUTION \n DESCRIPTION: Specify the verification solution OPTIONS: see \link Verification_Solution_Map \endlink \n DEFAULT: NO_VERIFICATION_SOLUTION \ingroup Config*/
  addEnumOption("KIND_VERIFICATION_SOLUTION", Kind_Verification_Solution, Verification_Solution_Map, VERIFICATION_SOLUTION::NONE);

  /*!\brief KIND_MATRIX_COLORING \n DESCRIPTION: Specify the method for matrix coloring for Jacobian computations OPTIONS: see \link MatrixColoring_Map \endlink \n DEFAULT GREEDY_COLORING \ingroup Config*/
  addEnumOption("KIND_MATRIX_COLORING", Kind_Matrix_Coloring, MatrixColoring_Map, GREEDY_COLORING);

  /*!\brief WEAKLY_COUPLED_HEAT_EQUATION \n DESCRIPTION: Enable heat equation for incompressible flows. \ingroup Config*/
  addBoolOption("WEAKLY_COUPLED_HEAT_EQUATION", Weakly_Coupled_Heat, NO);

  /*\brief AXISYMMETRIC \n DESCRIPTION: Axisymmetric simulation \n DEFAULT: false \ingroup Config */
  addBoolOption("AXISYMMETRIC", Axisymmetric, false);
  /* DESCRIPTION: Add the gravity force */
  addBoolOption("GRAVITY_FORCE", GravityForce, false);
  /* DESCRIPTION: Add the Vorticity Confinement term*/
  addBoolOption("VORTICITY_CONFINEMENT", VorticityConfinement, false);
  /* DESCRIPTION: Apply a body force as a source term (NO, YES) */
  addBoolOption("BODY_FORCE", Body_Force, false);
  body_force[0] = 0.0; body_force[1] = 0.0; body_force[2] = 0.0;
  /* DESCRIPTION: Vector of body force values (BodyForce_X, BodyForce_Y, BodyForce_Z) */
  addDoubleArrayOption("BODY_FORCE_VECTOR", 3, body_force);

  /* DESCRIPTION: Apply a body force as a source term for periodic boundary conditions \n Options: NONE, PRESSURE_DROP, MASSFLOW \n DEFAULT: NONE \ingroup Config */
  addEnumOption("KIND_STREAMWISE_PERIODIC", Kind_Streamwise_Periodic, Streamwise_Periodic_Map, ENUM_STREAMWISE_PERIODIC::NONE);
  /* DESCRIPTION: Use real periodicity for temperature \n Options: NO, YES \n DEFAULT: NO \ingroup Config */
  addBoolOption("STREAMWISE_PERIODIC_TEMPERATURE", Streamwise_Periodic_Temperature, false);
  /* DESCRIPTION: Heatflux boundary at streamwise periodic 'outlet', choose heat [W] such that net domain heatflux is zero. Only active if STREAMWISE_PERIODIC_TEMPERATURE is active. \n DEFAULT: 0.0 \ingroup Config */
  addDoubleOption("STREAMWISE_PERIODIC_OUTLET_HEAT", Streamwise_Periodic_OutletHeat, 0.0);
  /* DESCRIPTION: Delta pressure [Pa] on which basis body force will be computed, serves as initial value if MASSFLOW is chosen. \n DEFAULT: 1.0 \ingroup Config */
  addDoubleOption("STREAMWISE_PERIODIC_PRESSURE_DROP", Streamwise_Periodic_PressureDrop, 1.0);
  /* DESCRIPTION: Target Massflow [kg/s], Delta P will be adapted until m_dot is met. \n DEFAULT: 0.0 \ingroup Config  */
  addDoubleOption("STREAMWISE_PERIODIC_MASSFLOW", Streamwise_Periodic_TargetMassFlow, 0.0);

  /*!\brief RESTART_SOL \n DESCRIPTION: Restart solution from native solution file \n Options: NO, YES \ingroup Config */
  addBoolOption("RESTART_SOL", Restart, false);
  /*!\brief WRT_RESTART_COMPACT \n DESCRIPTION: Minimize the size of restart files \n Options: NO, YES \ingroup Config */
  addBoolOption("WRT_RESTART_COMPACT", Wrt_Restart_Compact, true);
  /*!\brief BINARY_RESTART \n DESCRIPTION: Read binary SU2 native restart files. \n Options: YES, NO \ingroup Config */
  addBoolOption("READ_BINARY_RESTART", Read_Binary_Restart, true);
  /*!\brief WRT_RESTART_OVERWRITE \n DESCRIPTION: overwrite restart files or append iteration number. \n Options: YES, NO \ingroup Config */
  addBoolOption("WRT_RESTART_OVERWRITE", Wrt_Restart_Overwrite, true);
  /*!\brief WRT_SURFACE_OVERWRITE \n DESCRIPTION: overwrite visualisation files or append iteration number. \n Options: YES, NO \ingroup Config */
  addBoolOption("WRT_SURFACE_OVERWRITE", Wrt_Surface_Overwrite, true);
  /*!\brief WRT_VOLUME_OVERWRITE \n DESCRIPTION: overwrite visualisation files or append iteration number. \n Options: YES, NO \ingroup Config */
  addBoolOption("WRT_VOLUME_OVERWRITE", Wrt_Volume_Overwrite, true);
  /*!\brief SYSTEM_MEASUREMENTS \n DESCRIPTION: System of measurements \n OPTIONS: see \link Measurements_Map \endlink \n DEFAULT: SI \ingroup Config*/
  addEnumOption("SYSTEM_MEASUREMENTS", SystemMeasurements, Measurements_Map, SI);

  /*!\par CONFIG_CATEGORY: FluidModel \ingroup Config*/
  /*!\brief FLUID_MODEL \n DESCRIPTION: Fluid model \n OPTIONS: See \link FluidModel_Map \endlink \n DEFAULT: STANDARD_AIR \ingroup Config*/
  addEnumOption("FLUID_MODEL", Kind_FluidModel, FluidModel_Map, STANDARD_AIR);
  /*!\brief FLUID_NAME \n DESCRIPTION: Fluid name \n OPTIONS: see coolprop homepage \n DEFAULT: nitrogen \ingroup Config*/
  addStringOption("FLUID_NAME", FluidName, string("nitrogen"));

  /*!\par CONFIG_CATEGORY: Data-driven fluid model parameters \ingroup Config*/
  /*!\brief INTERPOLATION_METHOD \n DESCRIPTION: Interpolation method used to determine the thermodynamic state of the fluid. \n OPTIONS: See \link DataDrivenMethod_Map \endlink DEFAULT: MLP \ingroup Config*/
  addEnumOption("INTERPOLATION_METHOD",Kind_DataDriven_Method, DataDrivenMethod_Map, ENUM_DATADRIVEN_METHOD::LUT);
  /*!\brief FILENAME_INTERPOLATOR \n DESCRIPTION: Input file for the interpolation method. \n \ingroup Config*/
  addStringListOption("FILENAMES_INTERPOLATOR", n_Datadriven_files, DataDriven_Method_FileNames);
  /*!\brief DATADRIVEN_NEWTON_RELAXATION \n DESCRIPTION: Relaxation factor for Newton solvers in data-driven fluid model. \n \ingroup Config*/
  addDoubleOption("DATADRIVEN_NEWTON_RELAXATION", DataDriven_Relaxation_Factor, 0.05);

  /*!\brief CONFINEMENT_PARAM \n DESCRIPTION: Input Confinement Parameter for Vorticity Confinement*/
  addDoubleOption("CONFINEMENT_PARAM", Confinement_Param, 0.0);

  /*!\par CONFIG_CATEGORY: Freestream Conditions \ingroup Config*/
  /*--- Options related to freestream specification ---*/

  /*!\brief GAS_CONSTANT \n DESCRIPTION: Specific gas constant (287.058 J/kg*K (air), only for compressible flows) \ingroup Config*/
  addDoubleOption("GAS_CONSTANT", Gas_Constant, 287.058);
  /*!\brief GAMMA_VALUE  \n DESCRIPTION: Ratio of specific heats (1.4 (air), only for compressible flows) \ingroup Config*/
  addDoubleOption("GAMMA_VALUE", Gamma, 1.4);
  /*!\brief THERMODYNAMIC_PRESSURE  \n DESCRIPTION: Thermodynamics(operating) Pressure (101325 Pa), only for incompressible flows) \ingroup Config*/
  addDoubleOption("THERMODYNAMIC_PRESSURE", Pressure_Thermodynamic, 101325.0);
  /*!\brief CP_VALUE  \n DESCRIPTION: Specific heat at constant pressure, Cp (1004.703 J/kg*K (air), constant density incompressible fluids only) \ingroup Config*/
  addDoubleListOption("SPECIFIC_HEAT_CP", nSpecific_Heat_Cp, Specific_Heat_Cp);
  /*!\brief THERMAL_EXPANSION_COEFF  \n DESCRIPTION: Thermal expansion coefficient (0.00347 K^-1 (air), used for Boussinesq approximation for liquids/non-ideal gases) \ingroup Config*/
  addDoubleOption("THERMAL_EXPANSION_COEFF", Thermal_Expansion_Coeff, 0.00347);
  /*!\brief MOLECULAR_WEIGHT \n DESCRIPTION: Molecular weight for an incompressible ideal gas (28.96 g/mol (air) default) \ingroup Config*/
  addDoubleListOption("MOLECULAR_WEIGHT", nMolecular_Weight, Molecular_Weight);

  ///* DESCRIPTION: Specify if Mutation++ library is used */
  /*--- Reading gas model as string or integer depending on TC library used. ---*/
  /* DESCRIPTION: Specify chemical model for multi-species simulations - read by Mutation++ library*/
  addStringOption("GAS_MODEL", GasModel, string("N2"));
  /* DESCRIPTION: Specify transport coefficient model for multi-species simulations */
  addEnumOption("TRANSPORT_COEFF_MODEL", Kind_TransCoeffModel, TransCoeffModel_Map, TRANSCOEFFMODEL::WILKE);
  /* DESCRIPTION: Specify mass fraction of each species */
  addDoubleListOption("GAS_COMPOSITION", nSpecies, Gas_Composition);
  /* DESCRIPTION: Specify mass fraction of each species for NEMO inlet*/
  addDoubleListOption("INLET_GAS_COMPOSITION", nSpecies_inlet, Inlet_MassFrac);
  /*!\brief INLET_TEMPERATURE_VE \n DESCRIPTION: NEMO inlet temperature_ve (K), if left 0 K, set to Ttr value \ingroup Config*/
  addDoubleOption("INLET_TEMPERATURE_VE", Inlet_Temperature_ve, 0.0);
  /* DESCRIPTION: Specify if mixture is frozen */
  addBoolOption("FROZEN_MIXTURE", frozen, false);
  /* DESCRIPTION: Specify if there is ionization */
  addBoolOption("IONIZATION", ionization, false);
  /* DESCRIPTION: Specify if there is VT transfer residual limiting */
  addBoolOption("VT_RESIDUAL_LIMITING", vt_transfer_res_limit, false);
  /* DESCRIPTION: List of catalytic walls */
  addStringListOption("CATALYTIC_WALL", nWall_Catalytic, Wall_Catalytic);
  /* DESCRIPTION: Specfify super-catalytic wall */
  addBoolOption("SUPERCATALYTIC_WALL", Supercatalytic_Wall, false);
  /* DESCRIPTION: Wall mass fractions for supercatalytic case */
  addDoubleListOption("SUPERCATALYTIC_WALL_COMPOSITION", nSpecies_Cat_Wall, Supercatalytic_Wall_Composition);
  /* DESCRIPTION: Specfify catalytic efficiency of wall if using gamma model */
  addDoubleOption("CATALYTIC_EFFICIENCY", CatalyticEfficiency, 1.0);
  /*!\brief MARKER_MONITORING\n DESCRIPTION: Marker(s) of the surface where evaluate the non-dimensional coefficients \ingroup Config*/

  /*--- Options related to VAN der WAALS MODEL and PENG ROBINSON ---*/

  /* DESCRIPTION: Critical Temperature, default value for AIR */
  addDoubleOption("CRITICAL_TEMPERATURE", Temperature_Critical, 131.00);
  /* DESCRIPTION: Critical Pressure, default value for MDM */
  addDoubleOption("CRITICAL_PRESSURE", Pressure_Critical, 3588550.0);
  /* DESCRIPTION: Critical Density, default value for MDM */
  addDoubleOption("CRITICAL_DENSITY", Density_Critical, 263.0);

  /*--- Options related to VAN der WAALS MODEL and PENG ROBINSON ---*/
  /* DESCRIPTION: Critical Density, default value for MDM */
   addDoubleOption("ACENTRIC_FACTOR", Acentric_Factor, 0.035);

   /*--- Options related to Viscosity Model ---*/
  /*!\brief VISCOSITY_MODEL \n DESCRIPTION: model of the viscosity \n OPTIONS: See \link ViscosityModel_Map \endlink \n DEFAULT: SUTHERLAND \ingroup Config*/
  addEnumOption("VISCOSITY_MODEL", Kind_ViscosityModel, ViscosityModel_Map, VISCOSITYMODEL::SUTHERLAND);

  /*--- Options related to Constant Viscosity Model ---*/

  /* DESCRIPTION: default value for AIR */
  addDoubleListOption("MU_CONSTANT", nMu_Constant, Mu_Constant);

  /*--- Options related to Sutherland Viscosity Model ---*/

  /* DESCRIPTION: Sutherland Viscosity Ref default value for AIR SI */
  addDoubleListOption("MU_REF", nMu_Ref, Mu_Ref);
  /* DESCRIPTION: Sutherland Temperature Ref, default value for AIR SI */
  addDoubleListOption("MU_T_REF", nMu_Temperature_Ref, Mu_Temperature_Ref);
  /* DESCRIPTION: Sutherland constant, default value for AIR SI */
  addDoubleListOption("SUTHERLAND_CONSTANT", nMu_S, Mu_S);

  /*--- Options related to Viscosity Model ---*/
  /*!\brief MIXINGVISCOSITY_MODEL \n DESCRIPTION: Mixing model of the viscosity \n OPTIONS: See \link ViscosityModel_Map \endlink \n DEFAULT: DAVIDSON \ingroup Config*/
  addEnumOption("MIXING_VISCOSITY_MODEL", Kind_MixingViscosityModel, MixingViscosityModel_Map, MIXINGVISCOSITYMODEL::DAVIDSON);

  /*--- Options related to Thermal Conductivity Model ---*/

  addEnumOption("CONDUCTIVITY_MODEL", Kind_ConductivityModel, ConductivityModel_Map, CONDUCTIVITYMODEL::CONSTANT_PRANDTL);

  /* DESCRIPTION: Definition of the turbulent thermal conductivity model (CONSTANT_PRANDTL_TURB (default), NONE). */
  addEnumOption("TURBULENT_CONDUCTIVITY_MODEL", Kind_ConductivityModel_Turb, TurbConductivityModel_Map, CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL);

 /*--- Options related to Constant Thermal Conductivity Model ---*/

 /* DESCRIPTION: default value for AIR */
  addDoubleListOption("THERMAL_CONDUCTIVITY_CONSTANT", nThermal_Conductivity_Constant , Thermal_Conductivity_Constant);

  /*--- Options related to temperature polynomial coefficients for fluid models. ---*/

  /* DESCRIPTION: Definition of the temperature polynomial coefficients for specific heat Cp. */
  addDoubleArrayOption("CP_POLYCOEFFS", N_POLY_COEFFS, cp_polycoeffs.data());
  /* DESCRIPTION: Definition of the temperature polynomial coefficients for specific heat Cp. */
  addDoubleArrayOption("MU_POLYCOEFFS", N_POLY_COEFFS, mu_polycoeffs.data());
  /* DESCRIPTION: Definition of the temperature polynomial coefficients for specific heat Cp. */
  addDoubleArrayOption("KT_POLYCOEFFS", N_POLY_COEFFS, kt_polycoeffs.data());

  /*!\brief REYNOLDS_NUMBER \n DESCRIPTION: Reynolds number (non-dimensional, based on the free-stream values). Needed for viscous solvers. For incompressible solvers the Reynolds length will always be 1.0 \n DEFAULT: 0.0 \ingroup Config */
  addDoubleOption("REYNOLDS_NUMBER", Reynolds, 0.0);
  /*!\brief REYNOLDS_LENGTH \n DESCRIPTION: Reynolds length (1 m by default). Used for compressible solver: incompressible solver will use 1.0. \ingroup Config */
  addDoubleOption("REYNOLDS_LENGTH", Length_Reynolds, 1.0);
  /*!\brief PRANDTL_LAM \n DESCRIPTION: Laminar Prandtl number (0.72 (air), only for compressible flows) \n DEFAULT: 0.72 \ingroup Config*/
  addDoubleListOption("PRANDTL_LAM", nPrandtl_Lam , Prandtl_Lam);
  /*!\brief PRANDTL_TURB \n DESCRIPTION: Turbulent Prandtl number (0.9 (air), only for compressible flows) \n DEFAULT 0.90 \ingroup Config*/
  addDoubleListOption("PRANDTL_TURB", nPrandtl_Turb , Prandtl_Turb);

  /*--- Options related to wall models. ---*/

  /*!\brief WALLMODEL_KAPPA \n DESCRIPTION: von Karman constant used for the wall model \n DEFAULT 0.41 \ingroup Config*/
  addDoubleOption("WALLMODEL_KAPPA", wallModel_Kappa, 0.41);
  /*!\brief WALLMODEL_MAXITER \n DESCRIPTION: Max iterations used for the wall model \n DEFAULT 200 \ingroup Config*/
  addUnsignedShortOption("WALLMODEL_MAXITER", wallModel_MaxIter, 200);
  /*!\brief WALLMODEL_RELFAC \n DESCRIPTION: Relaxation factor used for the wall model \n DEFAULT 0.5 \ingroup Config*/
  addDoubleOption("WALLMODEL_RELFAC", wallModel_RelFac, 0.5);
  /*!\brief WALLMODEL_MINYPLUS \n DESCRIPTION: lower limit for Y+ used for the wall model \n DEFAULT 5.0 \ingroup Config*/
  addDoubleOption("WALLMODEL_MINYPLUS", wallModel_MinYplus, 5.0);
  /*!\brief WALLMODEL_B \n DESCRIPTION: constant B used for the wall model \n DEFAULT 5.5 \ingroup Config*/
  addDoubleOption("WALLMODEL_B", wallModel_B, 5.5);

  /*!\brief BULK_MODULUS \n DESCRIPTION: Value of the Bulk Modulus  \n DEFAULT 1.42E5 \ingroup Config*/
  addDoubleOption("BULK_MODULUS", Bulk_Modulus, 1.42E5);
  /* DESCRIPTION: Epsilon^2 multipier in Beta calculation for incompressible preconditioner.  */
  addDoubleOption("BETA_FACTOR", Beta_Factor, 4.1);
  /*!\brief MACH_NUMBER  \n DESCRIPTION:  Mach number (non-dimensional, based on the free-stream values). 0.0 by default \ingroup Config*/
  addDoubleOption("MACH_NUMBER", Mach, 0.0);
  /*!\brief INIT_OPTION \n DESCRIPTION: Init option to choose between Reynolds or thermodynamics quantities for initializing the solution \n OPTIONS: see \link InitOption_Map \endlink \n DEFAULT REYNOLDS \ingroup Config*/
  addEnumOption("INIT_OPTION", Kind_InitOption, InitOption_Map, REYNOLDS);
  /* DESCRIPTION: Free-stream option to choose between density and temperature for initializing the solution */
  addEnumOption("FREESTREAM_OPTION", Kind_FreeStreamOption, FreeStreamOption_Map, FREESTREAM_OPTION::TEMPERATURE_FS);
  /*!\brief FREESTREAM_PRESSURE\n DESCRIPTION: Free-stream pressure (101325.0 N/m^2 by default) \ingroup Config*/
  addDoubleOption("FREESTREAM_PRESSURE", Pressure_FreeStream, 101325.0);
  /*!\brief FREESTREAM_DENSITY\n DESCRIPTION: Free-stream density (1.2886 Kg/m^3 (air), 998.2 Kg/m^3 (water)) \n DEFAULT -1.0 (calculated from others) \ingroup Config*/
  addDoubleOption("FREESTREAM_DENSITY", Density_FreeStream, -1.0);
  /*!\brief FREESTREAM_TEMPERATURE\n DESCRIPTION: Free-stream temperature (288.15 K by default) \ingroup Config*/
  addDoubleOption("FREESTREAM_TEMPERATURE", Temperature_FreeStream, 288.15);
  /*!\brief FREESTREAM_TEMPERATURE_VE\n DESCRIPTION: Free-stream vibrational-electronic temperature (288.15 K by default) \ingroup Config*/
  addDoubleOption("FREESTREAM_TEMPERATURE_VE", Temperature_ve_FreeStream, 288.15);


  /*--- Options related to incompressible flow solver ---*/

  /* DESCRIPTION: Option to choose the density model used in the incompressible flow solver. */
  addEnumOption("INC_DENSITY_MODEL", Kind_DensityModel, DensityModel_Map, INC_DENSITYMODEL::CONSTANT);
    /*!\brief ENERGY_EQUATION \n DESCRIPTION: Solve the energy equation in the incompressible flow solver. \ingroup Config*/
  addBoolOption("INC_ENERGY_EQUATION", Energy_Equation, false);
  /*!\brief INC_DENSITY_REF \n DESCRIPTION: Reference density for incompressible flows  \ingroup Config*/
  addDoubleOption("INC_DENSITY_REF", Inc_Density_Ref, 1.0);
  /*!\brief INC_VELOCITY_REF \n DESCRIPTION: Reference velocity for incompressible flows (1.0 by default) \ingroup Config*/
  addDoubleOption("INC_VELOCITY_REF", Inc_Velocity_Ref, 1.0);
  /*!\brief INC_TEMPERATURE_REF \n DESCRIPTION: Reference temperature for incompressible flows with the energy equation (1.0 by default) \ingroup Config*/
  addDoubleOption("INC_TEMPERATURE_REF", Inc_Temperature_Ref, 1.0);
  /*!\brief INC_DENSITY_INIT \n DESCRIPTION: Initial density for incompressible flows (1.2886 kg/m^3 by default) \ingroup Config*/
  addDoubleOption("INC_DENSITY_INIT", Inc_Density_Init, 1.2886);
  /*!\brief INC_VELOCITY_INIT \n DESCRIPTION: Initial velocity for incompressible flows (1.0,0,0 m/s by default) \ingroup Config*/
  vel_init[0] = 1.0; vel_init[1] = 0.0; vel_init[2] = 0.0;
  addDoubleArrayOption("INC_VELOCITY_INIT", 3, vel_init);
  /*!\brief INC_TEMPERATURE_INIT \n DESCRIPTION: Initial temperature for incompressible flows with the energy equation (288.15 K by default) \ingroup Config*/
  addDoubleOption("INC_TEMPERATURE_INIT", Inc_Temperature_Init, 288.15);
  /*!\brief INC_NONDIM \n DESCRIPTION: Non-dimensionalization scheme for incompressible flows. \ingroup Config*/
  addEnumOption("INC_NONDIM", Ref_Inc_NonDim, NonDim_Map, INITIAL_VALUES);
    /*!\brief INC_INLET_USENORMAL \n DESCRIPTION: Use the local boundary normal for the flow direction with the incompressible pressure inlet. \ingroup Config*/
  addBoolOption("INC_INLET_USENORMAL", Inc_Inlet_UseNormal, false);
  /*!\brief INC_INLET_DAMPING \n DESCRIPTION: Damping factor applied to the iterative updates to the velocity at a pressure inlet in incompressible flow (0.1 by default). \ingroup Config*/
  addDoubleOption("INC_INLET_DAMPING", Inc_Inlet_Damping, 0.1);
  /*!\brief INC_OUTLET_DAMPING \n DESCRIPTION: Damping factor applied to the iterative updates to the pressure at a mass flow outlet in incompressible flow (0.1 by default). \ingroup Config*/
  addDoubleOption("INC_OUTLET_DAMPING", Inc_Outlet_Damping, 0.1);

  /*--- Options related to the species solver. ---*/

  /*!\brief SPECIES_INIT \n DESCRIPTION: Initial values for scalar transport \ingroup Config*/
  addDoubleListOption("SPECIES_INIT", nSpecies_Init, Species_Init);
  /*!\brief SPECIES_CLIPPING \n DESCRIPTION: Activate clipping for scalar transport equations \n DEFAULT: false \ingroup Config*/
  addBoolOption("SPECIES_CLIPPING", Species_Clipping, false);
  /*!\brief SPECIES_CLIPPING_MAX \n DESCRIPTION: Maximum values for scalar clipping \ingroup Config*/
  addDoubleListOption("SPECIES_CLIPPING_MAX", nSpecies_Clipping_Max, Species_Clipping_Max);
  /*!\brief SPECIES_CLIPPING_MIN \n DESCRIPTION: Minimum values for scalar clipping \ingroup Config*/
  addDoubleListOption("SPECIES_CLIPPING_MIN", nSpecies_Clipping_Min, Species_Clipping_Min);

  /*!\brief FLAME_INIT_METHOD \n DESCRIPTION: Ignition method for flamelet solver \n DEFAULT: no ignition; cold flow only. */
  addEnumOption("FLAME_INIT_METHOD", flame_init_type, Flamelet_Init_Map, FLAMELET_INIT_TYPE::NONE);
  /*!\brief FLAME_INIT \n DESCRIPTION: flame front initialization using the flamelet model \ingroup Config*/
  /*--- flame offset (x,y,z) ---*/
  flame_init[0] = 0.0; flame_init[1] = 0.0; flame_init[2] = 0.0;
  /*--- flame normal (nx, ny, nz) ---*/
  flame_init[3] = 1.0; flame_init[4] = 0.0; flame_init[5] = 0.0;
  /*--- flame thickness (x) and flame burnt thickness (after this thickness, we have unburnt conditions again)  ---*/
  flame_init[6] = 0.5e-3; flame_init[7] = 1.0;
  addDoubleArrayOption("FLAME_INIT", 8,flame_init.begin());

  /*!\brief SPARK_INIT \n DESCRIPTION: spark initialization using the flamelet model \ingroup Config*/
  for (auto iSpark=0u; iSpark<6; ++iSpark) spark_init[iSpark]=0;
  addDoubleArrayOption("SPARK_INIT", 6, spark_init.begin());

  /*!\brief SPARK_REACTION_RATES \n DESCRIPTION: Net source term values applied to species within spark area during spark ignition. \ingroup Config*/
  addDoubleListOption("SPARK_REACTION_RATES", nspark, spark_reaction_rates);

  /*--- Options related to mass diffusivity and thereby the species solver. ---*/

  /*!\brief DIFFUSIVITY_MODEL\n DESCRIPTION: mass diffusivity model \n DEFAULT constant disffusivity \ingroup Config*/
  addEnumOption("DIFFUSIVITY_MODEL", Kind_Diffusivity_Model, Diffusivity_Model_Map, DIFFUSIVITYMODEL::CONSTANT_DIFFUSIVITY);
  /*!\brief DIFFUSIVITY_CONSTANT\n DESCRIPTION: mass diffusivity if DIFFUSIVITYMODEL::CONSTANT_DIFFUSIVITY is chosen \n DEFAULT 0.001 (Air) \ingroup Config*/
  addDoubleOption("DIFFUSIVITY_CONSTANT", Diffusivity_Constant , 0.001);
  /*!\brief SCHMIDT_LAM \n DESCRIPTION: Laminar Schmidt number of mass diffusion \n DEFAULT 1.0 (~for Gases) \ingroup Config*/
  addDoubleOption("SCHMIDT_NUMBER_LAMINAR", Schmidt_Number_Laminar, 1.0);
  /*!\brief SCHMIDT_TURB \n DESCRIPTION: Turbulent Schmidt number of mass diffusion \n DEFAULT 0.70 (more or less experimental value) \ingroup Config*/
  addDoubleOption("SCHMIDT_NUMBER_TURBULENT", Schmidt_Number_Turbulent, 0.7);
  /*!\brief DESCRIPTION: Constant Lewis number for mass diffusion */
  addDoubleListOption("CONSTANT_LEWIS_NUMBER", nConstant_Lewis_Number, Constant_Lewis_Number);

  vel_inf[0] = 1.0; vel_inf[1] = 0.0; vel_inf[2] = 0.0;
  /*!\brief FREESTREAM_VELOCITY\n DESCRIPTION: Free-stream velocity (m/s) */
  addDoubleArrayOption("FREESTREAM_VELOCITY", 3, vel_inf);
  /* DESCRIPTION: Free-stream viscosity (1.853E-5 Ns/m^2 (air), 0.798E-3 Ns/m^2 (water)) */
  addDoubleOption("FREESTREAM_VISCOSITY", Viscosity_FreeStream, -1.0);
  /* DESCRIPTION:  */
  addDoubleOption("FREESTREAM_INTERMITTENCY", Intermittency_FreeStream, 1.0);
  /* DESCRIPTION:  */
  addDoubleOption("FREESTREAM_TURBULENCEINTENSITY", TurbIntensityAndViscRatioFreeStream[0], 0.05);
  /* DESCRIPTION:  */
  addDoubleOption("N_CRITICAL", N_Critcal, 0.0);
  /* DESCRIPTION:  */
  addDoubleOption("FREESTREAM_NU_FACTOR", NuFactor_FreeStream, 3.0);
  /* DESCRIPTION:  */
  addDoubleOption("LOWER_LIMIT_K_FACTOR", KFactor_LowerLimit, 1.0e-15);
  /* DESCRIPTION:  */
  addDoubleOption("LOWER_LIMIT_OMEGA_FACTOR", OmegaFactor_LowerLimit, 1e-05);
  /* DESCRIPTION:  */
  addDoubleOption("ENGINE_NU_FACTOR", NuFactor_Engine, 3.0);
  /* DESCRIPTION:  */
  addDoubleOption("ACTDISK_SECONDARY_FLOW", SecondaryFlow_ActDisk, 0.0);
  /* DESCRIPTION:  */
  addDoubleOption("INITIAL_BCTHRUST", Initial_BCThrust, 4000.0);
  /* DESCRIPTION:  */
  addDoubleOption("FREESTREAM_TURB2LAMVISCRATIO", TurbIntensityAndViscRatioFreeStream[1], 10.0);
  /* DESCRIPTION: Side-slip angle (degrees, only for compressible flows) */
  addDoubleOption("SIDESLIP_ANGLE", AoS, 0.0);
  /*!\brief AOA  \n DESCRIPTION: Angle of attack (degrees, only for compressible flows) \ingroup Config*/
  addDoubleOption("AOA", AoA, 0.0);
  /* DESCRIPTION: Activate fixed CL mode (specify a CL instead of AoA). */
  addBoolOption("FIXED_CL_MODE", Fixed_CL_Mode, false);
  /* DESCRIPTION: Evaluate the dOF_dCL or dOF_dCMy during run time. */
  addBoolOption("EVAL_DOF_DCX", Eval_dOF_dCX, false);
  /* DESCRIPTION: DIscard the angle of attack in the solution and the increment in the geometry files. */
  addBoolOption("DISCARD_INFILES", Discard_InFiles, false);
  /* DESCRIPTION: Specify a fixed coefficient of lift instead of AoA (only for compressible flows) */
  addDoubleOption("TARGET_CL", Target_CL, 0.0);
  /* DESCRIPTION: Damping factor for fixed CL mode. */
  addDoubleOption("DCL_DALPHA", dCL_dAlpha, 0.2);
  /* DESCRIPTION: Maximum number of iterations between AoA updates for fixed CL problem. */
  addUnsignedLongOption("UPDATE_AOA_ITER_LIMIT", Update_AoA_Iter_Limit, 200);
  /* DESCRIPTION: Number of times Alpha is updated in a fix CL problem. */
  addUnsignedLongOption("UPDATE_IH", Update_iH, 5);
  /* DESCRIPTION: Number of iterations to evaluate dCL_dAlpha . */
  addUnsignedLongOption("ITER_DCL_DALPHA", Iter_dCL_dAlpha, 500);
  /* DESCRIPTION:  Value of dNetThrust/dBCThrust */
  addDoubleOption("DNETTHRUST_DBCTHRUST", dNetThrust_dBCThrust, 1.0);
  /* DESCRIPTION: Number of times Alpha is updated in a fix CL problem. */
  addUnsignedLongOption("UPDATE_BCTHRUST", Update_BCThrust, 5);


  /*!\par CONFIG_CATEGORY: Reference Conditions \ingroup Config*/
  /*--- Options related to reference values for nondimensionalization ---*/

  Length_Ref = 1.0; //<---- NOTE: this should be given an option or set as a const

  /*!\brief REF_ORIGIN_MOMENT_X\n DESCRIPTION: X Reference origin for moment computation \ingroup Config*/
  addDoubleListOption("REF_ORIGIN_MOMENT_X", nRefOriginMoment_X, RefOriginMoment_X);
  /*!\brief REF_ORIGIN_MOMENT_Y\n DESCRIPTION: Y Reference origin for moment computation \ingroup Config*/
  addDoubleListOption("REF_ORIGIN_MOMENT_Y", nRefOriginMoment_Y, RefOriginMoment_Y);
  /*!\brief REF_ORIGIN_MOMENT_Z\n DESCRIPTION: Z Reference origin for moment computation \ingroup Config*/
  addDoubleListOption("REF_ORIGIN_MOMENT_Z", nRefOriginMoment_Z, RefOriginMoment_Z);
  /*!\brief REF_AREA\n DESCRIPTION: Reference area for force coefficients (0 implies automatic calculation) \ingroup Config*/
  addDoubleOption("REF_AREA", RefArea, 1.0);
  /*!\brief SEMI_SPAN\n DESCRIPTION: Wing semi-span (0 implies automatic calculation) \ingroup Config*/
  addDoubleOption("SEMI_SPAN", SemiSpan, 0.0);
  /*!\brief REF_LENGTH\n DESCRIPTION: Reference length for pitching, rolling, and yawing non-dimensional moment \ingroup Config*/
  addDoubleOption("REF_LENGTH", RefLength, 1.0);
  /*!\brief REF_SHARP_EDGES\n DESCRIPTION: Reference coefficient for detecting sharp edges \ingroup Config*/
  addDoubleOption("REF_SHARP_EDGES", RefSharpEdges, 3.0);
  /*!\brief REF_VELOCITY\n DESCRIPTION: Reference velocity (incompressible only)  \ingroup Config*/
  addDoubleOption("REF_VELOCITY", Velocity_Ref, -1.0);
  /* !\brief REF_VISCOSITY  \n DESCRIPTION: Reference viscosity (incompressible only)  \ingroup Config*/
  addDoubleOption("REF_VISCOSITY", Viscosity_Ref, 1.0);
  /* DESCRIPTION: Type of mesh motion */
  addEnumOption("REF_DIMENSIONALIZATION", Ref_NonDim, NonDim_Map, DIMENSIONAL);

  /*!\par CONFIG_CATEGORY: Boundary Markers \ingroup Config*/
  /*--- Options related to various boundary markers ---*/

  /*!\brief MARKER_PLOTTING\n DESCRIPTION: Marker(s) of the surface in the surface flow solution file  \ingroup Config*/
  addStringListOption("MARKER_PLOTTING", nMarker_Plotting, Marker_Plotting);
  /*!\brief MARKER_MONITORING\n DESCRIPTION: Marker(s) of the surface where evaluate the non-dimensional coefficients \ingroup Config*/
  addStringListOption("MARKER_MONITORING", nMarker_Monitoring, Marker_Monitoring);

  /*!\brief MARKER_CONTROL_VOLUME\n DESCRIPTION: Marker(s) of the surface in the surface flow solution file  \ingroup Config*/
  addStringListOption("MARKER_ANALYZE", nMarker_Analyze, Marker_Analyze);
  /*!\brief MARKER_DESIGNING\n DESCRIPTION: Marker(s) of the surface where objective function (design problem) will be evaluated \ingroup Config*/
  addStringListOption("MARKER_DESIGNING", nMarker_Designing, Marker_Designing);
  /*!\brief GEO_MARKER\n DESCRIPTION: Marker(s) of the surface where evaluate the geometrical functions \ingroup Config*/
  addStringListOption("GEO_MARKER", nMarker_GeoEval, Marker_GeoEval);
  /*!\brief MARKER_EULER\n DESCRIPTION: Euler wall boundary marker(s) \ingroup Config*/
  addStringListOption("MARKER_EULER", nMarker_Euler, Marker_Euler);
  /*!\brief MARKER_FAR\n DESCRIPTION: Far-field boundary marker(s) \ingroup Config*/
  addStringListOption("MARKER_FAR", nMarker_FarField, Marker_FarField);
  /*!\brief MARKER_SYM\n DESCRIPTION: Symmetry boundary condition \ingroup Config*/
  addStringListOption("MARKER_SYM", nMarker_SymWall, Marker_SymWall);
  /*!\brief MARKER_NEARFIELD\n DESCRIPTION: Near-Field boundary condition \ingroup Config*/
  addStringListOption("MARKER_NEARFIELD", nMarker_NearFieldBound, Marker_NearFieldBound);
  /*!\brief MARKER_FLUID_INTERFACE\n DESCRIPTION: Fluid interface boundary marker(s) \ingroup Config*/
  addStringListOption("MARKER_FLUID_INTERFACE", nMarker_Fluid_InterfaceBound, Marker_Fluid_InterfaceBound);
  /*!\brief MARKER_DEFORM_MESH\n DESCRIPTION: Deformable marker(s) at the interface \ingroup Config*/
  addStringListOption("MARKER_DEFORM_MESH", nMarker_Deform_Mesh, Marker_Deform_Mesh);
  /*!\brief MARKER_DEFORM_MESH_SYM_PLANE\n DESCRIPTION: Symmetry plane for mesh deformation only \ingroup Config*/
  addStringListOption("MARKER_DEFORM_MESH_SYM_PLANE", nMarker_Deform_Mesh_Sym_Plane, Marker_Deform_Mesh_Sym_Plane);
  /*!\brief MARKER_FLUID_LOAD\n DESCRIPTION: Marker(s) in which the flow load is computed/applied \ingroup Config*/
  addStringListOption("MARKER_FLUID_LOAD", nMarker_Fluid_Load, Marker_Fluid_Load);
  /*!\brief MARKER_FSI_INTERFACE \n DESCRIPTION: ZONE interface boundary marker(s) \ingroup Config*/
  addStringListOption("MARKER_ZONE_INTERFACE", nMarker_ZoneInterface, Marker_ZoneInterface);
  /*!\brief MARKER_CHT_INTERFACE \n DESCRIPTION: CHT interface boundary marker(s) \ingroup Config*/
  addStringListOption("MARKER_CHT_INTERFACE", nMarker_CHTInterface, Marker_CHTInterface);
  /*!\brief CHT_INTERFACE_CONTACT_RESISTANCE: Thermal contact resistance values for each CHT inerface. \ingroup Config*/
  addDoubleListOption("CHT_INTERFACE_CONTACT_RESISTANCE", nMarker_ContactResistance, CHT_ContactResistance);
  /* DESCRIPTION: Internal boundary marker(s) */
  addStringListOption("MARKER_INTERNAL", nMarker_Internal, Marker_Internal);
  /* DESCRIPTION: Custom boundary marker(s) */
  addStringListOption("MARKER_CUSTOM", nMarker_Custom, Marker_Custom);
  /* DESCRIPTION: Periodic boundary marker(s)
   Format: ( periodic marker, donor marker, rotation_center_x, rotation_center_y,
   rotation_center_z, rotation_angle_x-axis, rotation_angle_y-axis,
   rotation_angle_z-axis, translation_x, translation_y, translation_z, ... ) */
  addPeriodicOption("MARKER_PERIODIC", nMarker_PerBound, Marker_PerBound, Marker_PerDonor,
                    Periodic_RotCenter, Periodic_RotAngles, Periodic_Translation);

  /*!\brief MARKER_PYTHON_CUSTOM\n DESCRIPTION: Python customizable marker(s) \ingroup Config*/
  addStringListOption("MARKER_PYTHON_CUSTOM", nMarker_PyCustom, Marker_PyCustom);

  /*!\brief MARKER_WALL_FUNCTIONS\n DESCRIPTION: Viscous wall markers for which wall functions must be applied.
   Format: (Wall function marker, wall function type, ...) \ingroup Config*/
  addWallFunctionOption("MARKER_WALL_FUNCTIONS", nMarker_WallFunctions, Marker_WallFunctions,
                        Kind_WallFunctions, IntInfo_WallFunctions, DoubleInfo_WallFunctions);

  /*!\brief MARKER_STRONG_BC\n DESCRIPTION: Markers where a strong BC must be applied.*/
  addStringListOption("MARKER_SPECIES_STRONG_BC", nMarker_StrongBC, Marker_StrongBC);

  /*!\brief ACTDISK_TYPE  \n DESCRIPTION: Actuator Disk boundary type \n OPTIONS: see \link ActDisk_Map \endlink \n Default: VARIABLES_JUMP \ingroup Config*/
  addEnumOption("ACTDISK_TYPE", Kind_ActDisk, ActDisk_Map, VARIABLES_JUMP);

  /*!\brief MARKER_ACTDISK\n DESCRIPTION: \ingroup Config*/
  addActDiskOption("MARKER_ACTDISK",
                   nMarker_ActDiskInlet, nMarker_ActDiskOutlet,  Marker_ActDiskInlet, Marker_ActDiskOutlet,
                   ActDisk_PressJump, ActDisk_TempJump, ActDisk_Omega);

  /*!\brief MARKER_ACTDISK_BEM_CG\n DESCRIPTION: Actuator disk CG for blade element momentum (BEM) method. \ingroup Config*/
  addActDiskBemOption("MARKER_ACTDISK_BEM_CG",
                      nMarker_ActDiskBemInlet_CG, nMarker_ActDiskBemOutlet_CG,  Marker_ActDiskBemInlet_CG, Marker_ActDiskBemOutlet_CG,
                      ActDiskBem_CG[0], ActDiskBem_CG[1], ActDiskBem_CG[2]);

  /*!\brief MARKER_ACTDISK_BEM_AXIS\n DESCRIPTION: Actuator disk axis for blade element momentum (BEM) method. \ingroup Config*/
  addActDiskBemOption("MARKER_ACTDISK_BEM_AXIS",
                      nMarker_ActDiskBemInlet_Axis, nMarker_ActDiskBemOutlet_Axis,  Marker_ActDiskBemInlet_Axis, Marker_ActDiskBemOutlet_Axis,
                      ActDiskBem_Axis[0], ActDiskBem_Axis[1], ActDiskBem_Axis[2]);

  /*!\brief ACTDISK_FILENAME \n DESCRIPTION: Input file for a specified actuator disk (w/ extension) \n DEFAULT: actdiskinput.dat \ingroup Config*/
  addStringOption("ACTDISK_FILENAME", ActDisk_FileName, string("actdiskinput.dat"));

  /*!\brief INLET_TYPE  \n DESCRIPTION: Inlet boundary type \n OPTIONS: see \link Inlet_Map \endlink \n DEFAULT: TOTAL_CONDITIONS \ingroup Config*/
  addEnumOption("INLET_TYPE", Kind_Inlet, Inlet_Map, INLET_TYPE::TOTAL_CONDITIONS);
  /*!\brief INC_INLET_TYPE \n DESCRIPTION: List of inlet types for incompressible flows. List length must match number of inlet markers. Options: VELOCITY_INLET, PRESSURE_INLET, INPUT_FILE. \ingroup Config*/
  addEnumListOption("INC_INLET_TYPE", nInc_Inlet, Kind_Inc_Inlet, Inlet_Map);
  addBoolOption("SPECIFIED_INLET_PROFILE", Inlet_From_File, false);
  /*!\brief INLET_FILENAME \n DESCRIPTION: Input file for a specified inlet profile (w/ extension) \n DEFAULT: inlet.dat \ingroup Config*/
  addStringOption("INLET_FILENAME", Inlet_Filename, string("inlet.dat"));
  /*!\brief INLET_MATCHING_TOLERANCE
   * \n DESCRIPTION: If a file is provided to specify the inlet profile,
   * this tolerance will be used to match the coordinates in the input file to
   * the points on the grid. \n DEFAULT: 1E-6 \ingroup Config*/
  addDoubleOption("INLET_MATCHING_TOLERANCE", Inlet_Matching_Tol, 1e-6);
  /*!\brief MARKER_INLET  \n DESCRIPTION: Inlet boundary marker(s) with the following formats,
   Total Conditions: (inlet marker, total temp, total pressure, flow_direction_x,
   flow_direction_y, flow_direction_z, ... ) where flow_direction is
   a unit vector.
   Mass Flow: (inlet marker, density, velocity magnitude, flow_direction_x,
   flow_direction_y, flow_direction_z, ... ) where flow_direction is
   a unit vector. \ingroup Config*/
  addInletOption("MARKER_INLET", nMarker_Inlet, Marker_Inlet, Inlet_Ttotal, Inlet_Ptotal, Inlet_FlowDir);
  /*!\brief MARKER_INLET_SPECIES \n DESCRIPTION: Inlet Species boundary marker(s) with the following format
   Inlet Species: (inlet_marker, Species1, Species2, ..., SpeciesN-1, inlet_marker2, Species1, Species2, ...) */
  addInletSpeciesOption("MARKER_INLET_SPECIES",nMarker_Inlet_Species, Marker_Inlet_Species, Inlet_SpeciesVal, nSpecies_per_Inlet);
  /*!\brief MARKER_INLET_TURBULENT \n DESCRIPTION: Inlet Turbulence boundary marker(s) with the following format
   Inlet Turbulent: (inlet_marker, TurbulentIntensity1, ViscosityRatio1, inlet_marker2, TurbulentIntensity2,
   ViscosityRatio2, ...) */
  addInletTurbOption("MARKER_INLET_TURBULENT", nMarker_Inlet_Turb, Marker_Inlet_Turb, Inlet_TurbVal, nTurb_Properties);
  /*!\brief MARKER_RIEMANN \n DESCRIPTION: Riemann boundary marker(s) with the following formats, a unit vector.
   * \n OPTIONS: See \link Riemann_Map \endlink. The variables indicated by the option and the flow direction unit vector must be specified. \ingroup Config*/
  addRiemannOption("MARKER_RIEMANN", nMarker_Riemann, Marker_Riemann, Kind_Data_Riemann, Riemann_Map, Riemann_Var1, Riemann_Var2, Riemann_FlowDir);
  /*!\brief MARKER_GILES \n DESCRIPTION: Giles boundary marker(s) with the following formats, a unit vector. */
  /* \n OPTIONS: See \link Giles_Map \endlink. The variables indicated by the option and the flow direction unit vector must be specified. \ingroup Config*/
  addGilesOption("MARKER_GILES", nMarker_Giles, Marker_Giles, Kind_Data_Giles, Giles_Map, Giles_Var1, Giles_Var2, Giles_FlowDir, RelaxFactorAverage, RelaxFactorFourier);
  /*!\brief SPATIAL_FOURIER \n DESCRIPTION: Option to compute the spatial fourier trasformation for the Giles BC. */
  addBoolOption("SPATIAL_FOURIER", SpatialFourier, false);
  /*!\brief GILES_EXTRA_RELAXFACTOR \n DESCRIPTION: the 1st coeff the value of the under relaxation factor to apply to the shroud and hub,
   * the 2nd coefficient is the the percentage of span-wise height influenced by this extra under relaxation factor.*/
  extrarelfac[0] = 0.1; extrarelfac[1] = 0.1;
  addDoubleArrayOption("GILES_EXTRA_RELAXFACTOR", 2, extrarelfac);
  /*!\brief AVERAGE_PROCESS_TYPE \n DESCRIPTION: types of mixing process for averaging quantities at the boundaries.
    \n OPTIONS: see \link MixingProcess_Map \endlink \n DEFAULT: AREA_AVERAGE \ingroup Config*/
  addEnumOption("MIXINGPLANE_INTERFACE_KIND", Kind_MixingPlaneInterface, MixingPlaneInterface_Map, NEAREST_SPAN);
  /*!\brief AVERAGE_PROCESS_KIND \n DESCRIPTION: types of mixing process for averaging quantities at the boundaries.
    \n OPTIONS: see \link MixingProcess_Map \endlink \n DEFAULT: AREA_AVERAGE \ingroup Config*/
  addEnumOption("AVERAGE_PROCESS_KIND", Kind_AverageProcess, AverageProcess_Map, AREA);
  /*!\brief PERFORMANCE_AVERAGE_PROCESS_KIND \n DESCRIPTION: types of mixing process for averaging quantities at the boundaries for performance computation.
      \n OPTIONS: see \link MixingProcess_Map \endlink \n DEFAULT: AREA_AVERAGE \ingroup Config*/
  addEnumOption("PERFORMANCE_AVERAGE_PROCESS_KIND", Kind_PerformanceAverageProcess, AverageProcess_Map, AREA);
  mixedout_coeff[0] = 1.0; mixedout_coeff[1] = 1.0E-05; mixedout_coeff[2] = 15.0;
  /*!\brief MIXEDOUT_COEFF \n DESCRIPTION: the 1st coeff is an under relaxation factor for the Newton method,
   * the 2nd coefficient is the tolerance for the Newton method, 3rd coefficient is the maximum number of
   * iteration for the Newton Method.*/
  addDoubleArrayOption("MIXEDOUT_COEFF", 3, mixedout_coeff);
  /*!\brief RAMP_ROTATING_FRAME\n DESCRIPTION: option to ramp up or down the rotating frame velocity value*/
  addBoolOption("RAMP_ROTATING_FRAME", RampRotatingFrame, false);
  rampRotFrame_coeff[0] = 0; rampRotFrame_coeff[1] = 1.0; rampRotFrame_coeff[2] = 1000.0;
      /*!\brief RAMP_ROTATING_FRAME_COEFF \n DESCRIPTION: the 1st coeff is the staring velocity,
   * the 2nd coeff is the number of iterations for the update, 3rd is the number of iteration */
  addDoubleArrayOption("RAMP_ROTATING_FRAME_COEFF", 3, rampRotFrame_coeff);
  /* DESCRIPTION: AVERAGE_MACH_LIMIT is a limit value for average procedure based on the mass flux. */
  addDoubleOption("AVERAGE_MACH_LIMIT", AverageMachLimit, 0.03);
  /*!\brief RAMP_OUTLET_PRESSURE\n DESCRIPTION: option to ramp up or down the rotating frame velocity value*/
  addBoolOption("RAMP_OUTLET_PRESSURE", RampOutletPressure, false);
  rampOutPres_coeff[0] = 100000.0; rampOutPres_coeff[1] = 1.0; rampOutPres_coeff[2] = 1000.0;
  /*!\brief RAMP_OUTLET_PRESSURE_COEFF \n DESCRIPTION: the 1st coeff is the staring outlet pressure,
   * the 2nd coeff is the number of iterations for the update, 3rd is the number of total iteration till reaching the final outlet pressure value */
  addDoubleArrayOption("RAMP_OUTLET_PRESSURE_COEFF", 3, rampOutPres_coeff);
  /*!\brief MARKER_MIXINGPLANE \n DESCRIPTION: Identify the boundaries in which the mixing plane is applied. \ingroup Config*/
  addStringListOption("MARKER_MIXINGPLANE_INTERFACE", nMarker_MixingPlaneInterface, Marker_MixingPlaneInterface);
  /*!\brief TURBULENT_MIXINGPLANE \n DESCRIPTION: Activate mixing plane also for turbulent quantities \ingroup Config*/
  addBoolOption("TURBULENT_MIXINGPLANE", turbMixingPlane, false);
  /*!\brief MARKER_TURBOMACHINERY \n DESCRIPTION: Identify the boundaries for which the turbomachinery settings are  applied. \ingroup Config*/
  addTurboPerfOption("MARKER_TURBOMACHINERY", nMarker_Turbomachinery, Marker_TurboBoundIn, Marker_TurboBoundOut, Marker_Turbomachinery);
  /*!\brief NUM_SPANWISE_SECTIONS \n DESCRIPTION: Integer number of spanwise sections to compute 3D turbo BC and Performance for turbomachinery */
  addUnsignedShortOption("NUM_SPANWISE_SECTIONS", nSpanWiseSections_User, 1);
  /*!\brief SPANWISE_KIND \n DESCRIPTION: type of algorithm to identify the span-wise sections at the turbo boundaries.
   \n OPTIONS: see \link SpanWise_Map \endlink \n Default: AUTOMATIC */
  addEnumOption("SPANWISE_KIND", Kind_SpanWise, SpanWise_Map, AUTOMATIC);
  /*!\brief TURBOMACHINERY_KIND \n DESCRIPTION: types of turbomachinery architecture.
      \n OPTIONS: see \link TurboMachinery_Map \endlink \n Default: AXIAL */
  addEnumListOption("TURBOMACHINERY_KIND",nTurboMachineryKind, Kind_TurboMachinery, TurboMachinery_Map);
  /*!\brief TURBOMACHINERY_KIND \n DESCRIPTION: types of turbomachynery Performance Calculations.
    \n OPTIONS: see \link TurboPerfKind_Map \endlink \n Default: TURBINE */
  addEnumListOption("TURBO_PERF_KIND", nTurboMachineryKind, Kind_TurboPerf, TurboPerfKind_Map);
  /*!\brief MARKER_SHROUD \n DESCRIPTION: markers in which velocity is forced to 0.0.
   * \n Format: (shroud1, shroud2, ...)*/
  addStringListOption("MARKER_SHROUD", nMarker_Shroud, Marker_Shroud);
  /*!\brief MARKER_SUPERSONIC_INLET  \n DESCRIPTION: Supersonic inlet boundary marker(s)
   * \n   Format: (inlet marker, temperature, static pressure, velocity_x,   velocity_y, velocity_z, ... ), i.e. primitive variables specified. \ingroup Config*/
  addInletOption("MARKER_SUPERSONIC_INLET", nMarker_Supersonic_Inlet, Marker_Supersonic_Inlet, Inlet_Temperature, Inlet_Pressure, Inlet_Velocity);
  /*!\brief MARKER_SUPERSONIC_OUTLET \n DESCRIPTION: Supersonic outlet boundary marker(s) \ingroup Config*/
  addStringListOption("MARKER_SUPERSONIC_OUTLET", nMarker_Supersonic_Outlet, Marker_Supersonic_Outlet);
  /*!\brief MARKER_OUTLET  \n DESCRIPTION: Outlet boundary marker(s)\n
   Format: ( outlet marker, back pressure (static), ... ) \ingroup Config*/
  addStringDoubleListOption("MARKER_OUTLET", nMarker_Outlet, Marker_Outlet, Outlet_Pressure);
  /*!\brief INC_INLET_TYPE \n DESCRIPTION: List of outlet types for incompressible flows. List length must match number of inlet markers. Options: PRESSURE_OUTLET, MASS_FLOW_OUTLET. \ingroup Config*/
  addEnumListOption("INC_OUTLET_TYPE", nInc_Outlet, Kind_Inc_Outlet, Inc_Outlet_Map);
  /*!\brief MARKER_ISOTHERMAL DESCRIPTION: Isothermal wall boundary marker(s)\n
   * Format: ( isothermal marker, wall temperature (static), ... ) \ingroup Config  */
  addStringDoubleListOption("MARKER_ISOTHERMAL", nMarker_Isothermal, Marker_Isothermal, Isothermal_Temperature);
  /*!\brief MARKER_HEATFLUX  \n DESCRIPTION: Specified heat flux wall boundary marker(s)
   Format: ( Heat flux marker, wall heat flux (static), ... ) \ingroup Config*/
  addStringDoubleListOption("MARKER_HEATFLUX", nMarker_HeatFlux, Marker_HeatFlux, Heat_Flux);
  /*!\brief INTEGRATED_HEATFLUX \n DESCRIPTION: Prescribe Heatflux in [W] instead of [W/m^2] \ingroup Config \default false */
  addBoolOption("INTEGRATED_HEATFLUX", Integrated_HeatFlux, false);
  /*!\brief MARKER_HEATTRANSFER DESCRIPTION: Heat flux with specified heat transfer coefficient boundary marker(s)\n
   * Format: ( Heat transfer marker, heat transfer coefficient, wall temperature (static), ... ) \ingroup Config  */
  addExhaustOption("MARKER_HEATTRANSFER", nMarker_HeatTransfer, Marker_HeatTransfer, HeatTransfer_Coeff, HeatTransfer_WallTemp);
  /*!\brief Smoluchowski/Maxwell wall boundary marker(s)  \n DESCRIPTION: Slip velocity and temperature jump wall boundary marker(s)
   Format: ( Heat flux marker,  wall temperature (static), momentum accomodation coefficient, thermal accomodation coefficient ... ) \ingroup Config*/
  addStringDoubleListOption("MARKER_SMOLUCHOWSKI_MAXWELL", nMarker_Smoluchowski_Maxwell, Marker_Smoluchowski_Maxwell, Isothermal_Temperature); //Missing TMAC and TAC
  /*!\brief WALL_ROUGHNESS  \n DESCRIPTION: Specified roughness heights at wall boundary marker(s)
   Format: ( Wall marker, roughness_height (static), ... ) \ingroup Config*/
  addStringDoubleListOption("WALL_ROUGHNESS", nRough_Wall, Marker_RoughWall, Roughness_Height);
  /*!\brief MARKER_ENGINE_INFLOW  \n DESCRIPTION: Engine inflow boundary marker(s)
   Format: ( nacelle inflow marker, fan face Mach, ... ) \ingroup Config*/
  addStringDoubleListOption("MARKER_ENGINE_INFLOW", nMarker_EngineInflow, Marker_EngineInflow, EngineInflow_Target);
  /* DESCRIPTION: Highlite area */
  addDoubleOption("HIGHLITE_AREA", Highlite_Area, 1.0);
  /* DESCRIPTION: Fan poly efficiency */
  addDoubleOption("FAN_POLY_EFF", Fan_Poly_Eff, 1.0);
  /*!\brief SUBSONIC_ENGINE\n DESCRIPTION: Engine subsonic intake region \ingroup Config*/
  addBoolOption("SUBSONIC_ENGINE", SubsonicEngine, false);
  /* DESCRIPTION: Actuator disk double surface */
  addBoolOption("ACTDISK_DOUBLE_SURFACE", ActDisk_DoubleSurface, false);

  /* DESCRIPTION: Actuator disk BEM switch for history file appending.*/
  addBoolOption("HISTORY_FILE_APPEND", History_File_Append_Flag, false);
  /* DESCRIPTION: Propeller blade angle for actuator disk BEM.*/
  addDoubleOption("BEM_PROP_BLADE_ANGLE", BEM_blade_angle, 23.9);
  /* DESCRIPTION: Propeller file name for actuator disk BEM.*/
  addStringOption("BEM_PROP_FILENAME", BEM_prop_filename, string("prop_geom_alfclcd_data.txt"));
  /* DESCRIPTION: Frequency for updating actuator disk with BEM.*/
  addUnsignedShortOption("BEM_FREQ", ActDiskBem_Frequency, 40);

  /* DESCRIPTION: Only half engine is in the computational grid */
  addBoolOption("ENGINE_HALF_MODEL", Engine_HalfModel, false);
  /* DESCRIPTION: Actuator disk SU2_DEF */
  addBoolOption("ACTDISK_SU2_DEF", ActDisk_SU2_DEF, false);
  /* DESCRIPTION: Definition of the distortion rack (radial number of proves / circumferential density (degree) */
  distortion[0] =  5.0; distortion[1] =  15.0;
  addDoubleArrayOption("DISTORTION_RACK", 2, distortion);
  /* DESCRIPTION: Values of the box to impose a subsonic nacellle (mach, Pressure, Temperature) */
  eng_val[0]=0.0; eng_val[1]=0.0; eng_val[2]=0.0; eng_val[3]=0.0;  eng_val[4]=0.0;
  addDoubleArrayOption("SUBSONIC_ENGINE_VALUES", 5, eng_val);
  /* DESCRIPTION: Coordinates of the box to impose a subsonic nacellle cylinder (Xmin, Ymin, Zmin, Xmax, Ymax, Zmax, Radius) */
  eng_cyl[0] = 0.0; eng_cyl[1] = 0.0; eng_cyl[2] = 0.0;
  eng_cyl[3] = 1E15; eng_cyl[4] = 1E15; eng_cyl[5] = 1E15; eng_cyl[6] = 1E15;
  addDoubleArrayOption("SUBSONIC_ENGINE_CYL", 7, eng_cyl);
  /* DESCRIPTION: Engine exhaust boundary marker(s)
   Format: (nacelle exhaust marker, total nozzle temp, total nozzle pressure, ... )*/
  addExhaustOption("MARKER_ENGINE_EXHAUST", nMarker_EngineExhaust, Marker_EngineExhaust, Exhaust_Temperature_Target, Exhaust_Pressure_Target);
  /* DESCRIPTION: Clamped boundary marker(s) */
  addStringListOption("MARKER_CLAMPED", nMarker_Clamped, Marker_Clamped);
  /* DESCRIPTION: Displacement boundary marker(s) */
  addStringDoubleListOption("MARKER_NORMAL_DISPL", nMarker_Displacement, Marker_Displacement, Displ_Value);
  /* DESCRIPTION: Load boundary marker(s) - uniform pressure in Pa */
  addStringDoubleListOption("MARKER_PRESSURE", nMarker_Load, Marker_Load, Load_Value);
  /* DESCRIPTION: Load boundary marker(s) */
  addStringDoubleListOption("MARKER_DAMPER", nMarker_Damper, Marker_Damper, Damper_Constant);
  /* DESCRIPTION: Load boundary marker(s)
   Format: (inlet marker, load, multiplier, dir_x, dir_y, dir_z, ... ), i.e. primitive variables specified. */
  addInletOption("MARKER_LOAD", nMarker_Load_Dir, Marker_Load_Dir, Load_Dir_Value, Load_Dir_Multiplier, Load_Dir);
  /* DESCRIPTION: Load boundary marker(s)
   Format: (inlet marker, load, multiplier, dir_x, dir_y, dir_z, ... ), i.e. primitive variables specified. */
  addInletOption("MARKER_DISPLACEMENT", nMarker_Disp_Dir, Marker_Disp_Dir, Disp_Dir_Value, Disp_Dir_Multiplier, Disp_Dir);
  /*!\brief SINE_LOAD\n DESCRIPTION: option to apply the load as a sine*/
  addBoolOption("SINE_LOAD", Sine_Load, false);
  sineload_coeff[0] = 0.0; sineload_coeff[1] = 0.0; sineload_coeff[2] = 0.0;
  /*!\brief SINE_LOAD_COEFF \n DESCRIPTION: the 1st coeff is the amplitude, the 2nd is the frequency, 3rd is the phase in radians */
  addDoubleArrayOption("SINE_LOAD_COEFF", 3, sineload_coeff);
  /*!\brief RAMP_AND_RELEASE\n DESCRIPTION: release the load after applying the ramp*/
  addBoolOption("RAMP_AND_RELEASE_LOAD", RampAndRelease, false);

  /* DESCRIPTION: Evaluation frequency for Engine and Actuator disk markers. */
  addUnsignedLongOption("BC_EVAL_FREQ", Bc_Eval_Freq, 40);
  /* DESCRIPTION: Damping factor for engine inlet condition */
  addDoubleOption("DAMP_ENGINE_INFLOW", Damp_Engine_Inflow, 0.95);
  /* DESCRIPTION: Damping factor for engine exhaust condition */
  addDoubleOption("DAMP_ENGINE_EXHAUST", Damp_Engine_Exhaust, 0.95);
  /*!\brief ENGINE_INFLOW_TYPE  \n DESCRIPTION: Inlet boundary type \n OPTIONS: see \link Engine_Inflow_Map \endlink \n Default: FAN_FACE_MACH \ingroup Config*/
  addEnumOption("ENGINE_INFLOW_TYPE", Kind_Engine_Inflow, Engine_Inflow_Map, FAN_FACE_MACH);
  /* DESCRIPTION: Evaluate a problem with engines */
  addBoolOption("ENGINE", Engine, false);

  /* DESCRIPTION:  Sharpness coefficient for the buffet sensor */
  addDoubleOption("BUFFET_K", Buffet_k, 10.0);
  /* DESCRIPTION:  Offset parameter for the buffet sensor */
  addDoubleOption("BUFFET_LAMBDA", Buffet_lambda, 0.0);

  /* DESCRIPTION: Use a Newton-Krylov method. */
  addBoolOption("NEWTON_KRYLOV", NewtonKrylov, false);
  /* DESCRIPTION: Integer parameters {startup iters, precond iters, initial tolerance relaxation}. */
  addUShortArrayOption("NEWTON_KRYLOV_IPARAM", NK_IntParam.size(), NK_IntParam.data());
  /* DESCRIPTION: Double parameters {startup residual drop, precond tolerance, full tolerance residual drop, findiff step}. */
  addDoubleArrayOption("NEWTON_KRYLOV_DPARAM", NK_DblParam.size(), NK_DblParam.data());

  /* DESCRIPTION: Number of samples for quasi-Newton methods. */
  addUnsignedShortOption("QUASI_NEWTON_NUM_SAMPLES", nQuasiNewtonSamples, 0);
  /* DESCRIPTION: Whether to use vectorized numerical schemes, less robust against transients. */
  addBoolOption("USE_VECTORIZATION", UseVectorization, false);

  /*!\par CONFIG_CATEGORY: Time-marching \ingroup Config*/
  /*--- Options related to time-marching ---*/

  /* DESCRIPTION: Unsteady simulation  */
  addEnumOption("TIME_MARCHING", TimeMarching, TimeMarching_Map, TIME_MARCHING::STEADY);
  /* DESCRIPTION:  Courant-Friedrichs-Lewy condition of the finest grid */
  addDoubleOption("CFL_NUMBER", CFLFineGrid, 1.25);
  /* DESCRIPTION:  Max time step in local time stepping simulations */
  addDoubleOption("MAX_DELTA_TIME", Max_DeltaTime, 1000000);
  /* DESCRIPTION: Activate The adaptive CFL number. */
  addBoolOption("CFL_ADAPT", CFL_Adapt, false);
  /* !\brief CFL_ADAPT_PARAM
   * DESCRIPTION: Parameters of the adaptive CFL number (factor down, factor up, CFL limit (min and max)[, acceptable linear residual][, starting iteration]).
   * Parameters in square brackets are optional, parameter "starting iteration" only valid with parameter "acceptable linear residual".
   * Factor down generally <1.0, factor up generally > 1.0 to cause the CFL to increase when the under-relaxation parameter is 1.0
   * and to decrease when the under-relaxation parameter is less than 0.1. Factor is multiplicative. \ingroup Config*/
  default_cfl_adapt[0] = 0.1;
  default_cfl_adapt[1] = 1.2;
  default_cfl_adapt[2] = 10.0;
  default_cfl_adapt[3] = 100.0;
  default_cfl_adapt[4] = 0.001;
  default_cfl_adapt[5] = 0.0;
  addDoubleListOption("CFL_ADAPT_PARAM", nCFL_AdaptParam, CFL_AdaptParam);
  /* DESCRIPTION: Reduction factor of the CFL coefficient in the adjoint problem */
  addDoubleOption("CFL_REDUCTION_ADJFLOW", CFLRedCoeff_AdjFlow, 0.8);
  /* DESCRIPTION: Reduction factor of the CFL coefficient in the level set problem */
  addDoubleOption("CFL_REDUCTION_TURB", CFLRedCoeff_Turb, 1.0);
  /* DESCRIPTION: Reduction factor of the CFL coefficient in the turbulent adjoint problem */
  addDoubleOption("CFL_REDUCTION_ADJTURB", CFLRedCoeff_AdjTurb, 1.0);
  /*!\brief CFL_REDUCTION_SPECIES \n DESCRIPTION: Reduction factor of the CFL coefficient in the species problem \n DEFAULT: 1.0 */
  addDoubleOption("CFL_REDUCTION_SPECIES", CFLRedCoeff_Species, 1.0);
  /* DESCRIPTION: External iteration offset due to restart */
  addUnsignedLongOption("EXT_ITER_OFFSET", ExtIter_OffSet, 0);
  // these options share nRKStep as their size, which is not a good idea in general
  /* DESCRIPTION: Runge-Kutta alpha coefficients */
  addDoubleListOption("RK_ALPHA_COEFF", nRKStep, RK_Alpha_Step);
  /* DESCRIPTION: Number of time levels for time accurate local time stepping. */
  addUnsignedShortOption("LEVELS_TIME_ACCURATE_LTS", nLevels_TimeAccurateLTS, 1);
  /* DESCRIPTION: Number of time DOFs used in the predictor step of ADER-DG. */
  addUnsignedShortOption("TIME_DOFS_ADER_DG", nTimeDOFsADER_DG, 2);
  /* DESCRIPTION: Unsteady Courant-Friedrichs-Lewy number of the finest grid */
  addDoubleOption("UNST_CFL_NUMBER", Unst_CFL, 0.0);
  /* DESCRIPTION: Integer number of periodic time instances for Harmonic Balance */
  addUnsignedShortOption("TIME_INSTANCES", nTimeInstances, 1);
  /* DESCRIPTION: Time period for Harmonic Balance wihtout moving meshes */
  addDoubleOption("HB_PERIOD", HarmonicBalance_Period, -1.0);
  /* DESCRIPTION:  Turn on/off harmonic balance preconditioning */
  addBoolOption("HB_PRECONDITION", HB_Precondition, false);
  /* DESCRIPTION: Starting direct solver iteration for the unsteady adjoint */
  addLongOption("UNST_ADJOINT_ITER", Unst_AdjointIter, 0);
  /* DESCRIPTION: Number of iterations to average the objective */
  addLongOption("ITER_AVERAGE_OBJ", Iter_Avg_Objective , 0);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_FLOW", Kind_TimeIntScheme_Flow, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_FEM_FLOW", Kind_TimeIntScheme_FEM_Flow, Time_Int_Map, RUNGE_KUTTA_EXPLICIT);
  /* DESCRIPTION: ADER-DG predictor step */
  addEnumOption("ADER_PREDICTOR", Kind_ADER_Predictor, Ader_Predictor_Map, ADER_ALIASED_PREDICTOR);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_ADJFLOW", Kind_TimeIntScheme_AdjFlow, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_TURB", Kind_TimeIntScheme_Turb, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_ADJTURB", Kind_TimeIntScheme_AdjTurb, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization for species equations */
  addEnumOption("TIME_DISCRE_SPECIES", Kind_TimeIntScheme_Species, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_FEA", Kind_TimeIntScheme_FEA, Time_Int_Map_FEA, STRUCT_TIME_INT::NEWMARK_IMPLICIT);
  /* DESCRIPTION: Time discretization for radiation problems*/
  addEnumOption("TIME_DISCRE_RADIATION", Kind_TimeIntScheme_Radiation, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIME_DISCRE_HEAT", Kind_TimeIntScheme_Heat, Time_Int_Map, EULER_IMPLICIT);
  /* DESCRIPTION: Time discretization */
  addEnumOption("TIMESTEP_HEAT", Kind_TimeStep_Heat, Heat_TimeStep_Map, MINIMUM);

  /*!\par CONFIG_CATEGORY: Linear solver definition \ingroup Config*/
  /*--- Options related to the linear solvers ---*/

  /*!\brief LINEAR_SOLVER
   *  \n DESCRIPTION: Linear solver for the implicit, mesh deformation, or discrete adjoint systems \n OPTIONS: see \link Linear_Solver_Map \endlink \n DEFAULT: FGMRES \ingroup Config*/
  addEnumOption("LINEAR_SOLVER", Kind_Linear_Solver, Linear_Solver_Map, FGMRES);
  /*!\brief LINEAR_SOLVER_PREC
   *  \n DESCRIPTION: Preconditioner for the Krylov linear solvers \n OPTIONS: see \link Linear_Solver_Prec_Map \endlink \n DEFAULT: LU_SGS \ingroup Config*/
  addEnumOption("LINEAR_SOLVER_PREC", Kind_Linear_Solver_Prec, Linear_Solver_Prec_Map, ILU);
  /* DESCRIPTION: Minimum error threshold for the linear solver for the implicit formulation */
  addDoubleOption("LINEAR_SOLVER_ERROR", Linear_Solver_Error, 1E-6);
  /* DESCRIPTION: Maximum number of iterations of the linear solver for the implicit formulation */
  addUnsignedLongOption("LINEAR_SOLVER_ITER", Linear_Solver_Iter, 10);
  /* DESCRIPTION: Fill in level for the ILU preconditioner */
  addUnsignedShortOption("LINEAR_SOLVER_ILU_FILL_IN", Linear_Solver_ILU_n, 0);
  /* DESCRIPTION: Maximum number of iterations of the linear solver for the implicit formulation */
  addUnsignedLongOption("LINEAR_SOLVER_RESTART_FREQUENCY", Linear_Solver_Restart_Frequency, 10);
  /* DESCRIPTION: Relaxation factor for iterative linear smoothers (SMOOTHER_ILU/JACOBI/LU-SGS/LINELET) */
  addDoubleOption("LINEAR_SOLVER_SMOOTHER_RELAXATION", Linear_Solver_Smoother_Relaxation, 1.0);
  /* DESCRIPTION: Custom number of threads used for additive domain decomposition for ILU and LU_SGS (0 is "auto"). */
  addUnsignedLongOption("LINEAR_SOLVER_PREC_THREADS", Linear_Solver_Prec_Threads, 0);
  /* DESCRIPTION: Relaxation factor for updates of adjoint variables. */
  addDoubleOption("RELAXATION_FACTOR_ADJOINT", Relaxation_Factor_Adjoint, 1.0);
  /* DESCRIPTION: Relaxation of the CHT coupling */
  addDoubleOption("RELAXATION_FACTOR_CHT", Relaxation_Factor_CHT, 1.0);
  /* DESCRIPTION: Roe coefficient */
  addDoubleOption("ROE_KAPPA", Roe_Kappa, 0.5);
  /* DESCRIPTION: Roe-Turkel preconditioning for low Mach number flows */
  addBoolOption("LOW_MACH_PREC", Low_Mach_Precon, false);
  /* DESCRIPTION: Post-reconstruction correction for low Mach number flows */
  addBoolOption("LOW_MACH_CORR", Low_Mach_Corr, false);
  /* DESCRIPTION: Minimum value for beta for the Roe-Turkel preconditioner */
  addDoubleOption("MIN_ROE_TURKEL_PREC", Min_Beta_RoeTurkel, 0.01);
  /* DESCRIPTION: Maximum value for beta for the Roe-Turkel preconditioner */
  addDoubleOption("MAX_ROE_TURKEL_PREC", Max_Beta_RoeTurkel, 0.2);
  /* DESCRIPTION: Entropy fix factor */
  addDoubleOption("ENTROPY_FIX_COEFF", EntropyFix_Coeff, 0.001);
  /* DESCRIPTION: Linear solver for the discete adjoint systems */
  addEnumOption("DISCADJ_LIN_SOLVER", Kind_DiscAdj_Linear_Solver, Linear_Solver_Map, FGMRES);
  /* DESCRIPTION: Preconditioner for the discrete adjoint Krylov linear solvers */
  addEnumOption("DISCADJ_LIN_PREC", Kind_DiscAdj_Linear_Prec, Linear_Solver_Prec_Map, ILU);
  /* DESCRIPTION: Linear solver for the discete adjoint systems */

  /*!\par CONFIG_CATEGORY: Convergence\ingroup Config*/
  /*--- Options related to convergence ---*/

  /*!\brief CONV_FIELD\n DESCRIPTION: Output field to monitor \n Default: depends on solver \ingroup Config*/
  addStringListOption("CONV_FIELD", nConvField, ConvField);
  /*!\brief CONV_RESIDUAL_MINVAL\n DESCRIPTION: Min value of the residual (log10 of the residual)\n DEFAULT: -14.0 \ingroup Config*/
  addDoubleOption("CONV_RESIDUAL_MINVAL", MinLogResidual, -14.0);
  /*!\brief CONV_STARTITER\n DESCRIPTION: Iteration number to begin convergence monitoring\n DEFAULT: 5 \ingroup Config*/
  addUnsignedLongOption("CONV_STARTITER", StartConv_Iter, 5);
  /*!\brief CONV_CAUCHY_ELEMS\n DESCRIPTION: Number of elements to apply the criteria. \n DEFAULT 100 \ingroup Config*/
  addUnsignedShortOption("CONV_CAUCHY_ELEMS", Cauchy_Elems, 100);
  /*!\brief CONV_CAUCHY_EPS\n DESCRIPTION: Epsilon to control the series convergence \n DEFAULT: 1e-10 \ingroup Config*/
  addDoubleOption("CONV_CAUCHY_EPS", Cauchy_Eps, 1E-10);

  /*!\brief CONV_WINDOW_STARTITER\n DESCRIPTION: Iteration number after START_ITER_WND  to begin convergence monitoring\n DEFAULT: 15 \ingroup Config*/
  addUnsignedLongOption("CONV_WINDOW_STARTITER", Wnd_StartConv_Iter, 15);
  /*!\brief CONV_WINDOW_CAUCHY_EPS\n DESCRIPTION: Epsilon to control the series convergence \n DEFAULT: 1e-3 \ingroup Config*/
  addDoubleOption("CONV_WINDOW_CAUCHY_EPS", Wnd_Cauchy_Eps, 1E-3);
  /*!\brief CONV_WINDOW_CAUCHY_ELEMS\n DESCRIPTION: Number of elements to apply the criteria. \n DEFAULT 100 \ingroup Config*/
  addUnsignedShortOption("CONV_WINDOW_CAUCHY_ELEMS", Wnd_Cauchy_Elems, 100);
  /*!\brief WINDOW_CAUCHY_CRIT \n DESCRIPTION: Determines, if the cauchy convergence criterion should be used for windowed time averaged objective functions*/
  addBoolOption("WINDOW_CAUCHY_CRIT",Wnd_Cauchy_Crit, false);
  /*!\brief CONV_WINDOW_FIELD
   * \n DESCRIPTION: Output fields  for the Cauchy criterium for the TIME iteration. The criterium is applied to the windowed time average of the chosen funcion. */
  addStringListOption("CONV_WINDOW_FIELD",nWndConvField, WndConvField);

  /*!\par CONFIG_CATEGORY: Multi-grid \ingroup Config*/
  /*!\brief MGLEVEL\n DESCRIPTION: Multi-grid Levels. DEFAULT: 0 \ingroup Config*/
  addUnsignedShortOption("MGLEVEL", nMGLevels, 0);
  /*!\brief MGCYCLE\n DESCRIPTION: Multi-grid cycle. OPTIONS: See \link MG_Cycle_Map \endlink. Defualt V_CYCLE \ingroup Config*/
  addEnumOption("MGCYCLE", MGCycle, MG_Cycle_Map, V_CYCLE);
  /*!\brief MG_PRE_SMOOTH\n DESCRIPTION: Multi-grid pre-smoothing level \ingroup Config*/
  addUShortListOption("MG_PRE_SMOOTH", nMG_PreSmooth, MG_PreSmooth);
  /*!\brief MG_POST_SMOOTH\n DESCRIPTION: Multi-grid post-smoothing level \ingroup Config*/
  addUShortListOption("MG_POST_SMOOTH", nMG_PostSmooth, MG_PostSmooth);
  /*!\brief MG_CORRECTION_SMOOTH\n DESCRIPTION: Jacobi implicit smoothing of the correction \ingroup Config*/
  addUShortListOption("MG_CORRECTION_SMOOTH", nMG_CorrecSmooth, MG_CorrecSmooth);
  /*!\brief MG_DAMP_RESTRICTION\n DESCRIPTION: Damping factor for the residual restriction. DEFAULT: 0.75 \ingroup Config*/
  addDoubleOption("MG_DAMP_RESTRICTION", Damp_Res_Restric, 0.75);
  /*!\brief MG_DAMP_PROLONGATION\n DESCRIPTION: Damping factor for the correction prolongation. DEFAULT 0.75 \ingroup Config*/
  addDoubleOption("MG_DAMP_PROLONGATION", Damp_Correc_Prolong, 0.75);

  /*!\par CONFIG_CATEGORY: Spatial Discretization \ingroup Config*/
  /*--- Options related to the spatial discretization ---*/

  /*!\brief NUM_METHOD_GRAD
   *  \n DESCRIPTION: Numerical method for spatial gradients \n OPTIONS: See \link Gradient_Map \endlink. \n DEFAULT: WEIGHTED_LEAST_SQUARES. \ingroup Config*/
  addEnumOption("NUM_METHOD_GRAD", Kind_Gradient_Method, Gradient_Map, WEIGHTED_LEAST_SQUARES);
  /*!\brief NUM_METHOD_GRAD
   *  \n DESCRIPTION: Numerical method for spatial gradients used only for upwind reconstruction \n OPTIONS: See \link Gradient_Map \endlink. \n DEFAULT: NO_GRADIENT. \ingroup Config*/
  addEnumOption("NUM_METHOD_GRAD_RECON", Kind_Gradient_Method_Recon, Gradient_Map, NO_GRADIENT);
  /*!\brief VENKAT_LIMITER_COEFF
   *  \n DESCRIPTION: Coefficient for the limiter. DEFAULT value 0.5. Larger values decrease the extent of limiting, values approaching zero cause lower-order approximation to the solution. \ingroup Config */
  addDoubleOption("VENKAT_LIMITER_COEFF", Venkat_LimiterCoeff, 0.05);
  /*!\brief ADJ_SHARP_LIMITER_COEFF
   *  \n DESCRIPTION: Coefficient for detecting the limit of the sharp edges. DEFAULT value 3.0.  Use with sharp edges limiter. \ingroup Config*/
  addDoubleOption("ADJ_SHARP_LIMITER_COEFF", AdjSharp_LimiterCoeff, 3.0);
  /*!\brief LIMITER_ITER
   *  \n DESCRIPTION: Freeze the value of the limiter after a number of iterations. DEFAULT value 999999. \ingroup Config*/
  addUnsignedLongOption("LIMITER_ITER", LimiterIter, 999999);

  /*!\brief CONV_NUM_METHOD_FLOW
   *  \n DESCRIPTION: Convective numerical method \n OPTIONS: See \link Upwind_Map \endlink , \link Centered_Map \endlink. \ingroup Config*/
  addConvectOption("CONV_NUM_METHOD_FLOW", Kind_ConvNumScheme_Flow, Kind_Centered_Flow, Kind_Upwind_Flow);

  /*!\brief NUM_METHOD_FEM_FLOW
   *  \n DESCRIPTION: Numerical method \n OPTIONS: See \link FEM_Map \endlink , \link Centered_Map \endlink. \ingroup Config*/
  addConvectFEMOption("NUM_METHOD_FEM_FLOW", Kind_ConvNumScheme_FEM_Flow, Kind_FEM_Flow);

  /*!\brief MUSCL_FLOW \n DESCRIPTION: Check if the MUSCL scheme should be used \ingroup Config*/
  addBoolOption("MUSCL_FLOW", MUSCL_Flow, true);
  /*!\brief SLOPE_LIMITER_FLOW
   * DESCRIPTION: Slope limiter for the direct solution. \n OPTIONS: See \link Limiter_Map \endlink \n DEFAULT VENKATAKRISHNAN \ingroup Config*/
  addEnumOption("SLOPE_LIMITER_FLOW", Kind_SlopeLimit_Flow, Limiter_Map, LIMITER::VENKATAKRISHNAN);
  jst_coeff[0] = 0.5; jst_coeff[1] = 0.02;
  /*!\brief JST_SENSOR_COEFF \n DESCRIPTION: 2nd and 4th order artificial dissipation coefficients for the JST method \ingroup Config*/
  addDoubleArrayOption("JST_SENSOR_COEFF", 2, jst_coeff);
  /*!\brief LAX_SENSOR_COEFF \n DESCRIPTION: 1st order artificial dissipation coefficients for the Lax-Friedrichs method. \ingroup Config*/
  addDoubleOption("LAX_SENSOR_COEFF", Kappa_1st_Flow, 0.15);
  /*!\brief USE_ACCURATE_FLUX_JACOBIANS \n DESCRIPTION: Use numerically computed Jacobians for AUSM+up(2) and SLAU(2) \ingroup Config*/
  addBoolOption("USE_ACCURATE_FLUX_JACOBIANS", Use_Accurate_Jacobians, false);
  /*!\brief CENTRAL_JACOBIAN_FIX_FACTOR \n DESCRIPTION: Improve the numerical properties (diagonal dominance) of the global Jacobian matrix, 3 to 4 is "optimum" (central schemes) \ingroup Config*/
  addDoubleOption("CENTRAL_JACOBIAN_FIX_FACTOR", Cent_Jac_Fix_Factor, 4.0);
  /*!\brief CENTRAL_JACOBIAN_FIX_FACTOR \n DESCRIPTION: Control numerical properties of the global Jacobian matrix using a multiplication factor for incompressible central schemes \ingroup Config*/
  addDoubleOption("CENTRAL_INC_JACOBIAN_FIX_FACTOR", Cent_Inc_Jac_Fix_Factor, 1.0);

  /*!\brief CONV_NUM_METHOD_ADJFLOW
   *  \n DESCRIPTION: Convective numerical method for the adjoint solver.
   *  \n OPTIONS:  See \link Upwind_Map \endlink , \link Centered_Map \endlink. Note: not all methods are guaranteed to be implemented for the adjoint solver. \ingroup Config */
  addConvectOption("CONV_NUM_METHOD_ADJFLOW", Kind_ConvNumScheme_AdjFlow, Kind_Centered_AdjFlow, Kind_Upwind_AdjFlow);
  /*!\brief MUSCL_ADJFLOW \n DESCRIPTION: Check if the MUSCL scheme should be used \ingroup Config*/
  addBoolOption("MUSCL_ADJFLOW", MUSCL_AdjFlow, true);
  /*!\brief SLOPE_LIMITER_ADJFLOW
     * DESCRIPTION: Slope limiter for the adjoint solution. \n OPTIONS: See \link Limiter_Map \endlink \n DEFAULT VENKATAKRISHNAN \ingroup Config*/
  addEnumOption("SLOPE_LIMITER_ADJFLOW", Kind_SlopeLimit_AdjFlow, Limiter_Map, LIMITER::VENKATAKRISHNAN);
  jst_adj_coeff[0] = 0.5; jst_adj_coeff[1] = 0.02;
  /*!\brief ADJ_JST_SENSOR_COEFF \n DESCRIPTION: 2nd and 4th order artificial dissipation coefficients for the adjoint JST method. \ingroup Config*/
  addDoubleArrayOption("ADJ_JST_SENSOR_COEFF", 2, jst_adj_coeff);
  /*!\brief LAX_SENSOR_COEFF \n DESCRIPTION: 1st order artificial dissipation coefficients for the adjoint Lax-Friedrichs method. \ingroup Config*/
  addDoubleOption("ADJ_LAX_SENSOR_COEFF", Kappa_1st_AdjFlow, 0.15);

  /*!\brief MUSCL_TURB \n DESCRIPTION: Check if the MUSCL scheme should be used \ingroup Config*/
  addBoolOption("MUSCL_TURB", MUSCL_Turb, false);
  /*!\brief SLOPE_LIMITER_TURB
   *  \n DESCRIPTION: Slope limiter  \n OPTIONS: See \link Limiter_Map \endlink \n DEFAULT VENKATAKRISHNAN \ingroup Config*/
  addEnumOption("SLOPE_LIMITER_TURB", Kind_SlopeLimit_Turb, Limiter_Map, LIMITER::VENKATAKRISHNAN);
  /*!\brief CONV_NUM_METHOD_TURB
   *  \n DESCRIPTION: Convective numerical method \ingroup Config*/
  addConvectOption("CONV_NUM_METHOD_TURB", Kind_ConvNumScheme_Turb, Kind_Centered_Turb, Kind_Upwind_Turb);

  /*!\brief MUSCL_ADJTURB \n DESCRIPTION: Check if the MUSCL scheme should be used \ingroup Config*/
  addBoolOption("MUSCL_ADJTURB", MUSCL_AdjTurb, false);
  /*!\brief SLOPE_LIMITER_ADJTURB
   *  \n DESCRIPTION: Slope limiter \n OPTIONS: See \link Limiter_Map \endlink \n DEFAULT VENKATAKRISHNAN \ingroup Config */
  addEnumOption("SLOPE_LIMITER_ADJTURB", Kind_SlopeLimit_AdjTurb, Limiter_Map, LIMITER::VENKATAKRISHNAN);
  /*!\brief CONV_NUM_METHOD_ADJTURB\n DESCRIPTION: Convective numerical method for the adjoint/turbulent problem \ingroup Config*/
  addConvectOption("CONV_NUM_METHOD_ADJTURB", Kind_ConvNumScheme_AdjTurb, Kind_Centered_AdjTurb, Kind_Upwind_AdjTurb);

  /*!\brief MUSCL_SPECIES \n DESCRIPTION: Check if the MUSCL scheme should be used \n DEFAULT false \ingroup Config*/
  addBoolOption("MUSCL_SPECIES", MUSCL_Species, false);
  /*!\brief SLOPE_LIMITER_SPECIES \n DESCRIPTION: Slope limiter \n OPTIONS: See \link Limiter_Map \endlink \n DEFAULT NONE \ingroup Config*/
  addEnumOption("SLOPE_LIMITER_SPECIES", Kind_SlopeLimit_Species, Limiter_Map, LIMITER::NONE);
  /*!\brief CONV_NUM_METHOD_SPECIES \n DESCRIPTION: Convective numerical method for species transport \ingroup Config*/
  addConvectOption("CONV_NUM_METHOD_SPECIES", Kind_ConvNumScheme_Species, Kind_Centered_Species, Kind_Upwind_Species);

  /*!\brief MUSCL_HEAT \n DESCRIPTION: Check if the MUSCL scheme should be used \ingroup Config*/
  addBoolOption("MUSCL_HEAT", MUSCL_Heat, false);
  /*!\brief SLOPE_LIMITER_HEAT \n DESCRIPTION: Slope limiter \n OPTIONS: See \link Limiter_Map \endlink \n DEFAULT NONE \ingroup Config*/
  addEnumOption("SLOPE_LIMITER_HEAT", Kind_SlopeLimit_Heat, Limiter_Map, LIMITER::NONE);
  /*!\brief CONV_NUM_METHOD_HEAT \n DESCRIPTION: Convective numerical method */
  addConvectOption("CONV_NUM_METHOD_HEAT", Kind_ConvNumScheme_Heat, Kind_Centered_Heat, Kind_Upwind_Heat);

  /*!\par CONFIG_CATEGORY: Adjoint and Gradient \ingroup Config*/
  /*--- Options related to the adjoint and gradient ---*/

  /*!\brief LIMIT_ADJFLOW \n DESCRIPTION: Limit value for the adjoint variable.\n DEFAULT: 1E6. \ingroup Config*/
  addDoubleOption("LIMIT_ADJFLOW", AdjointLimit, 1E6);
  /*!\brief MG_ADJFLOW\n DESCRIPTION: Multigrid with the adjoint problem. \n Defualt: YES \ingroup Config*/
  addBoolOption("MG_ADJFLOW", MG_AdjointFlow, true);

  /*!\brief OBJECTIVE_WEIGHT  \n DESCRIPTION: Adjoint problem boundary condition weights. Applies scaling factor to objective(s) \ingroup Config*/
  addDoubleListOption("OBJECTIVE_WEIGHT", nObjW, Weight_ObjFunc);
  /*!\brief OBJECTIVE_FUNCTION \n DESCRIPTION: Adjoint problem boundary condition \n OPTIONS: see \link Objective_Map \endlink \n DEFAULT: DRAG_COEFFICIENT \ingroup Config*/
  addEnumListOption("OBJECTIVE_FUNCTION", nObj, Kind_ObjFunc, Objective_Map);

  /*!\brief CUSTOM_OBJFUNC \n DESCRIPTION: User-provided definition of a custom objective function. \ingroup Config*/
  addStringOption("CUSTOM_OBJFUNC", CustomObjFunc, "");
  /*!\brief CUSTOM_OUTPUTS \n DESCRIPTION: User-provided definitions for custom output. \ingroup Config*/
  addStringOption("CUSTOM_OUTPUTS", CustomOutputs, "");

  /* DESCRIPTION: parameter for the definition of a complex objective function */
  addDoubleOption("DCD_DCL_VALUE", dCD_dCL, 0.0);
  /* DESCRIPTION: parameter for the definition of a complex objective function */
  addDoubleOption("DCMX_DCL_VALUE", dCMx_dCL, 0.0);
  /* DESCRIPTION: parameter for the definition of a complex objective function */
  addDoubleOption("DCMY_DCL_VALUE", dCMy_dCL, 0.0);
  /* DESCRIPTION: parameter for the definition of a complex objective function */
  addDoubleOption("DCMZ_DCL_VALUE", dCMz_dCL, 0.0);

  geo_loc[0] = 0.0; geo_loc[1] = 1.0;
  /* DESCRIPTION: Definition of the airfoil section */
  addDoubleArrayOption("GEO_BOUNDS", 2, geo_loc);
  /* DESCRIPTION: Identify the body to slice */
  addEnumOption("GEO_DESCRIPTION", Geo_Description, Geo_Description_Map, WING);
  /* DESCRIPTION: Z location of the waterline */
  addDoubleOption("GEO_WATERLINE_LOCATION", Geo_Waterline_Location, 0.0);
  /* DESCRIPTION: Number of section cuts to make when calculating internal volume */
  addUnsignedShortOption("GEO_NUMBER_STATIONS", nWingStations, 25);
  /* DESCRIPTION: Definition of the airfoil sections */
  addDoubleListOption("GEO_LOCATION_STATIONS", nLocationStations, LocationStations);
  nacelle_location[0] = 0.0; nacelle_location[1] = 0.0; nacelle_location[2] = 0.0;
  nacelle_location[3] = 0.0; nacelle_location[4] = 0.0;
  /* DESCRIPTION: Definition of the nacelle location (higlite coordinates, tilt angle, toe angle) */
  addDoubleArrayOption("GEO_NACELLE_LOCATION", 5, nacelle_location);
  /* DESCRIPTION: Output sectional forces for specified markers. */
  addBoolOption("GEO_PLOT_STATIONS", Plot_Section_Forces, false);
  /* DESCRIPTION: Mode of the GDC code (analysis, or gradient) */
  addEnumOption("GEO_MODE", GeometryMode, GeometryMode_Map, FUNCTION);

  /* DESCRIPTION: Drag weight in sonic boom Objective Function (from 0.0 to 1.0) */
  addDoubleOption("DRAG_IN_SONICBOOM", WeightCd, 0.0);
  /* DESCRIPTION: Sensitivity smoothing */
  addEnumOption("SENS_SMOOTHING", Kind_SensSmooth, Sens_Smoothing_Map, NO_SMOOTH);
  /* DESCRIPTION: Continuous Adjoint frozen viscosity */
  addBoolOption("FROZEN_VISC_CONT", Frozen_Visc_Cont, true);
  /* DESCRIPTION: Discrete Adjoint frozen viscosity */
  addBoolOption("FROZEN_VISC_DISC", Frozen_Visc_Disc, false);
  /* DESCRIPTION: Discrete Adjoint frozen limiter */
  addBoolOption("FROZEN_LIMITER_DISC", Frozen_Limiter_Disc, false);
  /* DESCRIPTION: Use an inconsistent (primal/dual) discrete adjoint formulation */
  addBoolOption("INCONSISTENT_DISC", Inconsistent_Disc, false);
   /* DESCRIPTION:  */
  addDoubleOption("FIX_AZIMUTHAL_LINE", FixAzimuthalLine, 90.0);
  /*!\brief SENS_REMOVE_SHARP
   * \n DESCRIPTION: Remove sharp edges from the sensitivity evaluation  \n Format: SENS_REMOVE_SHARP = YES \n DEFAULT: NO \ingroup Config*/
  addBoolOption("SENS_REMOVE_SHARP", Sens_Remove_Sharp, false);

  /* DESCRIPTION: Automatically reorient elements that seem flipped */
  addBoolOption("REORIENT_ELEMENTS",ReorientElements, true);

  /*!\par CONFIG_CATEGORY: Sobolev Gradient Solver Parameters \ingroup Config */
  /*--- Options related to the Sobolev smoothing solver ---*/

  /* DESCRIPTION: Switch to activate gradient smoothing */
  addBoolOption("SMOOTH_GRADIENT",SmoothGradient, false);
  /* DESCRIPTION: Epsilon of the identity term in the Laplace Beltrami Operator */
  addDoubleOption("SMOOTHING_EPSILON1",SmoothingEps1, 1.0);
  /* DESCRIPTION: Epsilon of the Laplace term in the Laplace Beltrami Operator */
  addDoubleOption("SMOOTHING_EPSILON2",SmoothingEps2, 1.0);
  /* DESCRIPTION: Switch to calculate for each dimension separately */
  addBoolOption("SEPARATE_DIMENSIONS", SmoothSepDim, false);
  /* DESCRIPTION: Switch to activate working on the design surfaces only */
  addBoolOption("SMOOTH_ON_SURFACE",SmoothOnSurface, false);
  /* DESCRIPTION: Switch to activate zero Dirichlet boundary for surface mode */
  addBoolOption("DIRICHLET_SURFACE_BOUNDARY", SmoothDirichletSurfaceBound, false);
  /* DESCRIPTION: Switch to activate the debbuging modus */
  addEnumOption("SOBOLEV_MODE", SmoothNumMode, Sobolev_Modus_Map, ENUM_SOBOLEV_MODUS::NONE);
  /*!\brief HESS_OBJFUNC_FILENAME
   *  \n DESCRIPTION: Output filename for the Sobolev Hessian approximation.  \ingroup Config*/
  addStringOption("HESS_OBJFUNC_FILENAME", ObjFunc_Hess_FileName, string("of_hess.dat"));

  /*  DESCRIPTION: Linear solver for the gradient smoothing\n OPTIONS: see \link Linear_Solver_Map \endlink \n DEFAULT: FGMRES \ingroup Config*/
  addEnumOption("GRAD_LINEAR_SOLVER", Kind_Grad_Linear_Solver, Linear_Solver_Map, FGMRES);
  /*  \n DESCRIPTION: Preconditioner for the Krylov linear solvers \n OPTIONS: see \link Linear_Solver_Prec_Map \endlink \n DEFAULT: ILU \ingroup Config*/
  addEnumOption("GRAD_LINEAR_SOLVER_PREC", Kind_Grad_Linear_Solver_Prec, Linear_Solver_Prec_Map, ILU);
  /* DESCRIPTION: Minimum error threshold for the linear solver for the implicit formulation */
  addDoubleOption("GRAD_LINEAR_SOLVER_ERROR", Grad_Linear_Solver_Error, 1E-14);
  /* DESCRIPTION: Maximum number of iterations of the linear solver for the implicit formulation */
  addUnsignedLongOption("GRAD_LINEAR_SOLVER_ITER", Grad_Linear_Solver_Iter, 1000);

  /*!\par CONFIG_CATEGORY: Input/output files and formats \ingroup Config */
  /*--- Options related to input/output files and formats ---*/

  /*!\brief OUTPUT_FORMAT \n DESCRIPTION: I/O format for output plots. \n OPTIONS: see \link TabOutput_Map \endlink \n DEFAULT: TECPLOT \ingroup Config */
  addEnumOption("TABULAR_FORMAT", Tab_FileFormat, TabOutput_Map, TAB_OUTPUT::TAB_CSV);
  /*!\brief OUTPUT_PRECISION \n DESCRIPTION: Set <ofstream>.precision(value) to specified value for SU2_DOT and HISTORY output. Useful for exact gradient validation. \n DEFAULT: 6 \ingroup Config */
  addUnsignedShortOption("OUTPUT_PRECISION", output_precision, 10);
  /*!\brief ACTDISK_JUMP \n DESCRIPTION: The jump is given by the difference in values or a ratio */
  addEnumOption("ACTDISK_JUMP", ActDisk_Jump, Jump_Map, DIFFERENCE);
  /*!\brief MESH_FORMAT \n DESCRIPTION: Mesh input file format \n OPTIONS: see \link Input_Map \endlink \n DEFAULT: SU2 \ingroup Config*/
  addEnumOption("MESH_FORMAT", Mesh_FileFormat, Input_Map, SU2);
  /* DESCRIPTION:  Mesh input file */
  addStringOption("MESH_FILENAME", Mesh_FileName, string("mesh.su2"));
  /*!\brief MESH_OUT_FILENAME \n DESCRIPTION: Mesh output file name. Used when converting, scaling, or deforming a mesh. \n DEFAULT: mesh_out.su2 \ingroup Config*/
  addStringOption("MESH_OUT_FILENAME", Mesh_Out_FileName, string("mesh_out.su2"));

  /* DESCRIPTION: List of the number of grid points in the RECTANGLE or BOX grid in the x,y,z directions. (default: (33,33,33) ). */
  addShortListOption("MESH_BOX_SIZE", nMesh_Box_Size, Mesh_Box_Size);

  /* DESCRIPTION: List of the length of the RECTANGLE or BOX grid in the x,y,z directions. (default: (1.0,1.0,1.0) ).  */
  mesh_box_length[0] = 1.0; mesh_box_length[1] = 1.0; mesh_box_length[2] = 1.0;
  addDoubleArrayOption("MESH_BOX_LENGTH", 3, mesh_box_length);

  /* DESCRIPTION: List of the offset from 0.0 of the RECTANGLE or BOX grid in the x,y,z directions. (default: (0.0,0.0,0.0) ). */
  mesh_box_offset[0] = 0.0; mesh_box_offset[1] = 0.0; mesh_box_offset[2] = 0.0;
  addDoubleArrayOption("MESH_BOX_OFFSET", 3, mesh_box_offset);

  /* DESCRIPTION: Polynomial degree of the FEM solution for the RECTANGLE or BOX grid. (default: 1). */
  addUnsignedShortOption("MESH_BOX_POLY_SOL_FEM", Mesh_Box_PSolFEM, 1);

  /* DESCRIPTION: Determine if the mesh file supports multizone. \n DEFAULT: true (temporarily) */
  addBoolOption("MULTIZONE_MESH", Multizone_Mesh, true);
  /* DESCRIPTION: Determine if we need to allocate memory to store the multizone residual. \n DEFAULT: false (temporarily) */
  addBoolOption("MULTIZONE_RESIDUAL", Multizone_Residual, false);

  /* !\brief CONTROLLING_VARIABLE_NAMES \n DESCRIPTION: Names of the variables used as inputs for the data regression method in flamelet or data-driven fluid models. */
  addStringListOption("CONTROLLING_VARIABLE_NAMES", n_control_vars, controlling_variable_names);

  /* !\brief CONTROLLING_VARIABLE_SOURCE_NAMES \n DESCRIPTION: Names of the variables in the flamelet manifold corresponding to the source terms of the controlling variables. */
  addStringListOption("CONTROLLING_VARIABLE_SOURCE_NAMES", n_control_vars, cv_source_names);

  /* DESCRIPTION: Names of the passive lookup variables for flamelet LUT */
  addStringListOption("LOOKUP_NAMES", n_lookups, lookup_names);

  /* DESCRIPTION: Names of the user transport equations solved in the flamelet problem. */
  addStringListOption("USER_SCALAR_NAMES", n_user_scalars, user_scalar_names);

  /* DESCRIPTION: Names of the user scalar source terms. */
  addStringListOption("USER_SOURCE_NAMES", n_user_sources, user_source_names);

  /* DESCRIPTION: Enable preferential diffusion for FGM simulations. \n DEFAULT: false */
  addBoolOption("PREFERENTIAL_DIFFUSION", preferential_diffusion, false);

  /*!\brief CONV_FILENAME \n DESCRIPTION: Output file convergence history (w/o extension) \n DEFAULT: history \ingroup Config*/
  addStringOption("CONV_FILENAME", Conv_FileName, string("history"));
  /*!\brief BREAKDOWN_FILENAME \n DESCRIPTION: Output file forces breakdown \ingroup Config*/
  addStringOption("BREAKDOWN_FILENAME", Breakdown_FileName, string("forces_breakdown.dat"));
  /*!\brief SOLUTION_FLOW_FILENAME \n DESCRIPTION: Restart flow input file (the file output under the filename set by RESTART_FLOW_FILENAME) \n DEFAULT: solution_flow.dat \ingroup Config */
  addStringOption("SOLUTION_FILENAME", Solution_FileName, string("solution.dat"));
  /*!\brief SOLUTION_ADJ_FILENAME\n DESCRIPTION: Restart adjoint input file. Objective function abbreviation is expected. \ingroup Config*/
  addStringOption("SOLUTION_ADJ_FILENAME", Solution_AdjFileName, string("solution_adj.dat"));
  /*!\brief RESTART_FLOW_FILENAME \n DESCRIPTION: Output file restart flow \ingroup Config*/
  addStringOption("RESTART_FILENAME", Restart_FileName, string("restart.dat"));
  /*!\brief RESTART_ADJ_FILENAME  \n DESCRIPTION: Output file restart adjoint. Objective function abbreviation will be appended. \ingroup Config*/
  addStringOption("RESTART_ADJ_FILENAME", Restart_AdjFileName, string("restart_adj.dat"));
  /*!\brief VOLUME_FLOW_FILENAME  \n DESCRIPTION: Output file flow (w/o extension) variables \ingroup Config */
  addStringOption("VOLUME_FILENAME", Volume_FileName, string("vol_solution"));
  /*!\brief VOLUME_ADJ_FILENAME
   *  \n DESCRIPTION: Output file adjoint (w/o extension) variables  \ingroup Config*/
  addStringOption("VOLUME_ADJ_FILENAME", Adj_FileName, string("adj_vol_solution"));
  /*!\brief GRAD_OBJFUNC_FILENAME
   *  \n DESCRIPTION: Output objective function gradient  \ingroup Config*/
  addStringOption("GRAD_OBJFUNC_FILENAME", ObjFunc_Grad_FileName, string("of_grad.dat"));
  /*!\brief VALUE_OBJFUNC_FILENAME
   *  \n DESCRIPTION: Output objective function  \ingroup Config*/
  addStringOption("VALUE_OBJFUNC_FILENAME", ObjFunc_Value_FileName, string("of_func.dat"));
  /*!\brief SURFACE_FLOW_FILENAME
   *  \n DESCRIPTION: Output file surface flow coefficient (w/o extension)  \ingroup Config*/
  addStringOption("SURFACE_FILENAME", SurfCoeff_FileName, string("surface"));
  /*!\brief SURFACE_ADJ_FILENAME
   *  \n DESCRIPTION: Output file surface adjoint coefficient (w/o extension)  \ingroup Config*/
  addStringOption("SURFACE_ADJ_FILENAME", SurfAdjCoeff_FileName, string("surface_adjoint"));
  /*!\brief SURFACE_SENS_FILENAME_FILENAME
   *  \n DESCRIPTION: Output file surface sensitivity (discrete adjoint) (w/o extension)  \ingroup Config*/
  addStringOption("SURFACE_SENS_FILENAME", SurfSens_FileName, string("surface_sens"));
  /*!\brief VOLUME_SENS_FILENAME
   *  \n DESCRIPTION: Output file volume sensitivity (discrete adjoint))  \ingroup Config*/
  addStringOption("VOLUME_SENS_FILENAME", VolSens_FileName, string("volume_sens"));
  /* DESCRIPTION: Output the performance summary to the console at the end of SU2_CFD  \ingroup Config*/
  addBoolOption("WRT_PERFORMANCE", Wrt_Performance, false);
  /* DESCRIPTION: Output the tape statistics (discrete adjoint)  \ingroup Config*/
  addBoolOption("WRT_AD_STATISTICS", Wrt_AD_Statistics, false);
  /*!\brief MARKER_ANALYZE_AVERAGE
   *  \n DESCRIPTION: Output averaged flow values on specified analyze marker.
   *  Options: AREA, MASSFLUX
   *  \n Use with MARKER_ANALYZE. \ingroup Config*/
  addEnumOption("MARKER_ANALYZE_AVERAGE", Kind_Average, Average_Map, AVERAGE_MASSFLUX);
  /*!\brief COMM_LEVEL
   *  \n DESCRIPTION: Level of MPI communications during runtime  \ingroup Config*/
  addEnumOption("COMM_LEVEL", Comm_Level, Comm_Map, COMM_FULL);

  /*!\par CONFIG_CATEGORY: Dynamic mesh definition \ingroup Config*/
  /*--- Options related to dynamic meshes ---*/

  /* DESCRIPTION: Type of mesh motion */
  addEnumOption("GRID_MOVEMENT", Kind_GridMovement, GridMovement_Map, NO_MOVEMENT);
  /* DESCRIPTION: Type of surface motion */
  addEnumListOption("SURFACE_MOVEMENT",nKind_SurfaceMovement, Kind_SurfaceMovement, SurfaceMovement_Map);
  /* DESCRIPTION: Marker(s) of moving surfaces (MOVING_WALL or DEFORMING grid motion). */
  addStringListOption("MARKER_MOVING", nMarker_Moving, Marker_Moving);
  /* DESCRIPTION: Marker(s) of gradient problem boundaries. */
  addStringListOption("MARKER_SOBOLEVBC", nMarker_SobolevBC, Marker_SobolevBC);
  /* DESCRIPTION: Mach number (non-dimensional, based on the mesh velocity and freestream vals.) */
  addDoubleOption("MACH_MOTION", Mach_Motion, 0.0);
  /* DESCRIPTION: Coordinates of the rigid motion origin */
  addDoubleArrayOption("MOTION_ORIGIN", 3, Motion_Origin);
  /* DESCRIPTION: Translational velocity vector (m/s) in the x, y, & z directions (RIGID_MOTION only) */
  addDoubleArrayOption("TRANSLATION_RATE", 3, Translation_Rate);
  /* DESCRIPTION: Angular velocity vector (rad/s) about x, y, & z axes (RIGID_MOTION only) */
  addDoubleArrayOption("ROTATION_RATE", 3, Rotation_Rate);
  /* DESCRIPTION: Pitching angular freq. (rad/s) about x, y, & z axes (RIGID_MOTION only) */
  addDoubleArrayOption("PITCHING_OMEGA", 3, Pitching_Omega);
  /* DESCRIPTION: Pitching amplitude (degrees) about x, y, & z axes (RIGID_MOTION only) */
  addDoubleArrayOption("PITCHING_AMPL", 3, Pitching_Ampl);
  /* DESCRIPTION: Pitching phase offset (degrees) about x, y, & z axes (RIGID_MOTION only) */
  addDoubleArrayOption("PITCHING_PHASE", 3, Pitching_Phase);
  /* DESCRIPTION: Plunging angular freq. (rad/s) in x, y, & z directions (RIGID_MOTION only) */
  addDoubleArrayOption("PLUNGING_OMEGA", 3, Plunging_Omega);
  /* DESCRIPTION: Plunging amplitude (m) in x, y, & z directions (RIGID_MOTION only) */
  addDoubleArrayOption("PLUNGING_AMPL", 3, Plunging_Ampl);
  /* DESCRIPTION: Coordinates of the rigid motion origin */
  addDoubleListOption("SURFACE_MOTION_ORIGIN", nMarkerMotion_Origin, MarkerMotion_Origin);
  /* DESCRIPTION: Translational velocity vector (m/s) in the x, y, & z directions (DEFORMING only) */
  addDoubleListOption("SURFACE_TRANSLATION_RATE", nMarkerTranslation, MarkerTranslation_Rate);
  /* DESCRIPTION: Angular velocity vector (rad/s) about x, y, & z axes (DEFORMING only) */
  addDoubleListOption("SURFACE_ROTATION_RATE", nMarkerRotation_Rate, MarkerRotation_Rate);
  /* DESCRIPTION: Pitching angular freq. (rad/s) about x, y, & z axes (DEFORMING only) */
  addDoubleListOption("SURFACE_PITCHING_OMEGA", nMarkerPitching_Omega, MarkerPitching_Omega);
  /* DESCRIPTION: Pitching amplitude (degrees) about x, y, & z axes (DEFORMING only) */
  addDoubleListOption("SURFACE_PITCHING_AMPL", nMarkerPitching_Ampl, MarkerPitching_Ampl);
  /* DESCRIPTION: Pitching phase offset (degrees) about x, y, & z axes (DEFORMING only) */
  addDoubleListOption("SURFACE_PITCHING_PHASE", nMarkerPitching_Phase, MarkerPitching_Phase);
  /* DESCRIPTION: Plunging angular freq. (rad/s) in x, y, & z directions (DEFORMING only) */
  addDoubleListOption("SURFACE_PLUNGING_OMEGA", nMarkerPlunging_Omega, MarkerPlunging_Omega);
  /* DESCRIPTION: Plunging amplitude (m) in x, y, & z directions (DEFORMING only) */
  addDoubleListOption("SURFACE_PLUNGING_AMPL", nMarkerPlunging_Ampl, MarkerPlunging_Ampl);
  /* DESCRIPTION: Value to move motion origins (1 or 0) */
  addUShortListOption("MOVE_MOTION_ORIGIN", nMoveMotion_Origin, MoveMotion_Origin);

  /*!\par CONFIG_CATEGORY: Aeroelastic Simulation (Typical Section Model) \ingroup Config*/
  /*--- Options related to aeroelastic simulations using the Typical Section Model) ---*/
  /* DESCRIPTION: The flutter speed index (modifies the freestream condition) */
  addDoubleOption("FLUTTER_SPEED_INDEX", FlutterSpeedIndex, 0.6);
  /* DESCRIPTION: Natural frequency of the spring in the plunging direction (rad/s). */
  addDoubleOption("PLUNGE_NATURAL_FREQUENCY", PlungeNaturalFrequency, 100);
  /* DESCRIPTION: Natural frequency of the spring in the pitching direction (rad/s). */
  addDoubleOption("PITCH_NATURAL_FREQUENCY", PitchNaturalFrequency, 100);
  /* DESCRIPTION: The airfoil mass ratio. */
  addDoubleOption("AIRFOIL_MASS_RATIO", AirfoilMassRatio, 60);
  /* DESCRIPTION: Distance in semichords by which the center of gravity lies behind the elastic axis. */
  addDoubleOption("CG_LOCATION", CG_Location, 1.8);
  /* DESCRIPTION: The radius of gyration squared (expressed in semichords) of the typical section about the elastic axis. */
  addDoubleOption("RADIUS_GYRATION_SQUARED", RadiusGyrationSquared, 3.48);
  /* DESCRIPTION: Solve the aeroelastic equations every given number of internal iterations. */
  addUnsignedShortOption("AEROELASTIC_ITER", AeroelasticIter, 3);

  /*!\par CONFIG_CATEGORY: Optimization Problem*/

  /* DESCRIPTION: Scale the line search in the optimizer */
  addDoubleOption("OPT_RELAX_FACTOR", Opt_RelaxFactor, 1.0);

  /* DESCRIPTION: Bound the line search in the optimizer */
  addDoubleOption("OPT_LINE_SEARCH_BOUND", Opt_LineSearch_Bound, 1E6);

  /*!\par CONFIG_CATEGORY: Wind Gust \ingroup Config*/
  /*--- Options related to wind gust simulations ---*/

  /* DESCRIPTION: Apply a wind gust */
  addBoolOption("WIND_GUST", Wind_Gust, false);
  /* DESCRIPTION: Type of gust */
  addEnumOption("GUST_TYPE", Gust_Type, Gust_Type_Map, NO_GUST);
  /* DESCRIPTION: Gust wavelenght (meters) */
  addDoubleOption("GUST_WAVELENGTH", Gust_WaveLength, 0.0);
  /* DESCRIPTION: Number of gust periods */
  addDoubleOption("GUST_PERIODS", Gust_Periods, 1.0);
  /* DESCRIPTION: Gust amplitude (m/s) */
  addDoubleOption("GUST_AMPL", Gust_Ampl, 0.0);
  /* DESCRIPTION: Time at which to begin the gust (sec) */
  addDoubleOption("GUST_BEGIN_TIME", Gust_Begin_Time, 0.0);
  /* DESCRIPTION: Location at which the gust begins (meters) */
  addDoubleOption("GUST_BEGIN_LOC", Gust_Begin_Loc, 0.0);
  /* DESCRIPTION: Direction of the gust X or Y dir */
  addEnumOption("GUST_DIR", Gust_Dir, Gust_Dir_Map, Y_DIR);

  /* Fixed values for turbulence quantities to keep them at inflow conditions. */
  /* DESCRIPTION: Fix turbulence quantities to far-field values inside an upstream half-space. */
  addBoolOption("TURB_FIXED_VALUES", Turb_Fixed_Values, false);
  /* DESCRIPTION: Shift of the fixed values half-space, in length units in the direction of far-field velocity. */
  addDoubleOption("TURB_FIXED_VALUES_DOMAIN", Turb_Fixed_Values_MaxScalarProd, numeric_limits<su2double>::lowest());

  /* Harmonic Balance config */
  /* DESCRIPTION: Omega_HB = 2*PI*frequency - frequencies for Harmonic Balance method */
  addDoubleListOption("OMEGA_HB", nOmega_HB, Omega_HB);

  /*!\par CONFIG_CATEGORY: Equivalent Area \ingroup Config*/
  /*--- Options related to the equivalent area ---*/

  /* DESCRIPTION: Evaluate equivalent area on the Near-Field  */
  addBoolOption("EQUIV_AREA", EquivArea, false);
  ea_lim[0] = 0.0; ea_lim[1] = 1.0; ea_lim[2] = 1.0;
  /* DESCRIPTION: Integration limits of the equivalent area ( xmin, xmax, Dist_NearField ) */
  addDoubleArrayOption("EA_INT_LIMIT", 3, ea_lim);
  /* DESCRIPTION: Equivalent area scaling factor */
  addDoubleOption("EA_SCALE_FACTOR", EA_ScaleFactor, 1.0);

  // these options share nDV as their size in the option references; not a good idea
  /*!\par CONFIG_CATEGORY: Grid deformation \ingroup Config*/
  /*--- Options related to the grid deformation ---*/

  /* DESCRIPTION: Kind of deformation */
  addEnumListOption("DV_KIND", nDV, Design_Variable, Param_Map);
  /* DESCRIPTION: Marker of the surface to which we are going apply the shape deformation */
  addStringListOption("DV_MARKER", nMarker_DV, Marker_DV);
  /* DESCRIPTION: Parameters of the shape deformation
   - FFD_CONTROL_POINT_2D ( FFDBox ID, i_Ind, j_Ind, x_Disp, y_Disp )
   - FFD_RADIUS_2D ( FFDBox ID )
   - FFD_CAMBER_2D ( FFDBox ID, i_Ind )
   - FFD_THICKNESS_2D ( FFDBox ID, i_Ind )
   - HICKS_HENNE ( Lower Surface (0)/Upper Surface (1)/Only one Surface (2), x_Loc )
   - SURFACE_BUMP ( x_start, x_end, x_Loc )
   - CST ( Lower Surface (0)/Upper Surface (1), Kulfan parameter number, Total number of Kulfan parameters for surface )
   - NACA_4DIGITS ( 1st digit, 2nd digit, 3rd and 4th digit )
   - PARABOLIC ( Center, Thickness )
   - TRANSLATION ( x_Disp, y_Disp, z_Disp )
   - ROTATION ( x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
   - OBSTACLE ( Center, Bump size )
   - SPHERICAL ( ControlPoint_Index, Theta_Disp, R_Disp )
   - FFD_CONTROL_POINT ( FFDBox ID, i_Ind, j_Ind, k_Ind, x_Disp, y_Disp, z_Disp )
   - FFD_TWIST ( FFDBox ID, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
   - FFD_ROTATION ( FFDBox ID, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
   - FFD_CONTROL_SURFACE ( FFDBox ID, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
   - FFD_CAMBER ( FFDBox ID, i_Ind, j_Ind )
   - FFD_THICKNESS ( FFDBox ID, i_Ind, j_Ind ) */
  addDVParamOption("DV_PARAM", nDV, ParamDV, FFDTag, Design_Variable);
  /* DESCRIPTION: New value of the shape deformation */
  addDVValueOption("DV_VALUE", nDV_Value, DV_Value, nDV, ParamDV, Design_Variable);
  /* DESCRIPTION: Provide a file of surface positions from an external parameterization. */
  addStringOption("DV_FILENAME", DV_Filename, string("surface_positions.dat"));
  /* DESCRIPTION: File of sensitivities as an unordered ASCII file with rows of x, y, z, dJ/dx, dJ/dy, dJ/dz for each volume grid point. */
  addStringOption("DV_UNORDERED_SENS_FILENAME", DV_Unordered_Sens_Filename, string("unordered_sensitivity.dat"));
  /* DESCRIPTION: File of sensitivities as an ASCII file with rows of x, y, z, dJ/dx, dJ/dy, dJ/dz for each surface grid point. */
  addStringOption("DV_SENS_FILENAME", DV_Sens_Filename, string("surface_sensitivity.dat"));
  /*!\brief OUTPUT_FORMAT \n DESCRIPTION: I/O format for output plots. \n OPTIONS: see \link Output_Map \endlink \n DEFAULT: TECPLOT \ingroup Config */
  addEnumOption("DV_SENSITIVITY_FORMAT", Sensitivity_FileFormat, Sensitivity_Map, SU2_NATIVE);
  /* DESCRIPTION: Hold the grid fixed in a region */
  addBoolOption("HOLD_GRID_FIXED", Hold_GridFixed, false);
  grid_fix[0] = -1E15; grid_fix[1] = -1E15; grid_fix[2] = -1E15;
  grid_fix[3] =  1E15; grid_fix[4] =  1E15; grid_fix[5] =  1E15;
  /* DESCRIPTION: Coordinates of the box where the grid will be deformed (Xmin, Ymin, Zmin, Xmax, Ymax, Zmax) */
  addDoubleArrayOption("HOLD_GRID_FIXED_COORD", 6, grid_fix);

  /*!\par CONFIG_CATEGORY: Deformable mesh \ingroup Config*/
  /*--- option related to deformable meshes ---*/
  /* DESCRIPTION: Decide whether the mesh will undergo deformations */
  addBoolOption("DEFORM_MESH", Deform_Mesh, false);
  /* DESCRIPTION: Print the residuals during mesh deformation to the console */
  addBoolOption("DEFORM_CONSOLE_OUTPUT", Deform_Output, false);
  /* DESCRIPTION: Number of nonlinear deformation iterations (surface deformation increments) */
  addUnsignedLongOption("DEFORM_NONLINEAR_ITER", GridDef_Nonlinear_Iter, 1);
  /* DESCRIPTION: Deform coefficient (-1.0 to 0.5) */
  addDoubleOption("DEFORM_COEFF", Deform_Coeff, 1E6);
  /* DESCRIPTION: Deform limit in m or inches */
  addDoubleOption("DEFORM_LIMIT", Deform_Limit, 1E6);
  /* DESCRIPTION: Type of element stiffness imposed for FEA mesh deformation (INVERSE_VOLUME, WALL_DISTANCE, CONSTANT_STIFFNESS) */
  addEnumOption("DEFORM_STIFFNESS_TYPE", Deform_StiffnessType, Deform_Stiffness_Map, SOLID_WALL_DISTANCE);
  /* DESCRIPTION: Young's modulus for constant stiffness FEA method of grid deformation */
  addDoubleOption("DEFORM_ELASTICITY_MODULUS", Deform_ElasticityMod, 2E11);
  /* DESCRIPTION: Poisson's ratio for constant stiffness FEA method of grid deformation */
  addDoubleOption("DEFORM_POISSONS_RATIO", Deform_PoissonRatio, 0.3);
  /* DESCRIPTION: Size of the layer of highest stiffness for wall distance-based mesh stiffness */
  addDoubleOption("DEFORM_STIFF_LAYER_SIZE", Deform_StiffLayerSize, 0.0);
  /*  DESCRIPTION: Linear solver for the mesh deformation\n OPTIONS: see \link Linear_Solver_Map \endlink \n DEFAULT: FGMRES \ingroup Config*/
  addEnumOption("DEFORM_LINEAR_SOLVER", Kind_Deform_Linear_Solver, Linear_Solver_Map, FGMRES);
  /*  \n DESCRIPTION: Preconditioner for the Krylov linear solvers \n OPTIONS: see \link Linear_Solver_Prec_Map \endlink \n DEFAULT: LU_SGS \ingroup Config*/
  addEnumOption("DEFORM_LINEAR_SOLVER_PREC", Kind_Deform_Linear_Solver_Prec, Linear_Solver_Prec_Map, ILU);
  /* DESCRIPTION: Minimum error threshold for the linear solver for the implicit formulation */
  addDoubleOption("DEFORM_LINEAR_SOLVER_ERROR", Deform_Linear_Solver_Error, 1E-14);
  /* DESCRIPTION: Maximum number of iterations of the linear solver for the implicit formulation */
  addUnsignedLongOption("DEFORM_LINEAR_SOLVER_ITER", Deform_Linear_Solver_Iter, 1000);

  /*!\par CONFIG_CATEGORY: FEM flow solver definition \ingroup Config*/
  /*--- Options related to the finite element flow solver---*/

  /* DESCRIPTION: Riemann solver used for DG (ROE, LAX-FRIEDRICH, AUSM, HLLC, VAN_LEER) */
  addEnumOption("RIEMANN_SOLVER_FEM", Riemann_Solver_FEM, Upwind_Map, UPWIND::ROE);
  /* DESCRIPTION: Constant factor applied for quadrature with straight elements (2.0 by default) */
  addDoubleOption("QUADRATURE_FACTOR_STRAIGHT_FEM", Quadrature_Factor_Straight, 2.0);
  /* DESCRIPTION: Constant factor applied for quadrature with curved elements (3.0 by default) */
  addDoubleOption("QUADRATURE_FACTOR_CURVED_FEM", Quadrature_Factor_Curved, 3.0);
  /* DESCRIPTION: Factor applied during quadrature in time for ADER-DG. (2.0 by default) */
  addDoubleOption("QUADRATURE_FACTOR_TIME_ADER_DG", Quadrature_Factor_Time_ADER_DG, 2.0);
  /* DESCRIPTION: Factor for the symmetrizing terms in the DG FEM discretization (1.0 by default) */
  addDoubleOption("THETA_INTERIOR_PENALTY_DG_FEM", Theta_Interior_Penalty_DGFEM, 1.0);
  /* DESCRIPTION: Compute the entropy in the fluid model (YES, NO) */
  addBoolOption("COMPUTE_ENTROPY_FLUID_MODEL", Compute_Entropy, true);
  /* DESCRIPTION: Use the lumped mass matrix for steady DGFEM computations */
  addBoolOption("USE_LUMPED_MASSMATRIX_DGFEM", Use_Lumped_MassMatrix_DGFEM, false);
  /* DESCRIPTION: Only compute the exact Jacobian of the spatial discretization (NO, YES) */
  addBoolOption("JACOBIAN_SPATIAL_DISCRETIZATION_ONLY", Jacobian_Spatial_Discretization_Only, false);

  /* DESCRIPTION: Number of aligned bytes for the matrix multiplications. Multiple of 64. (128 by default) */
  addUnsignedShortOption("ALIGNED_BYTES_MATMUL", byteAlignmentMatMul, 128);

  /*!\par CONFIG_CATEGORY: FEA solver \ingroup Config*/
  /*--- Options related to the FEA solver ---*/

  /*!\brief FEA_FILENAME \n DESCRIPTION: Filename to input for element-based properties \n Default: default_element_properties.dat \ingroup Config */
  addStringOption("FEA_FILENAME", FEA_FileName, string("default_element_properties.dat"));
  /* DESCRIPTION: Determine if advanced features are used from the element-based FEA analysis (NO, YES = experimental) */
  addBoolOption("FEA_ADVANCED_MODE", FEAAdvancedMode, false);

  /* DESCRIPTION: Modulus of elasticity */
  addDoubleListOption("ELASTICITY_MODULUS", nElasticityMod, ElasticityMod);
  /* DESCRIPTION: Poisson ratio */
  addDoubleListOption("POISSON_RATIO", nPoissonRatio, PoissonRatio);
  /* DESCRIPTION: Material density */
  addDoubleListOption("MATERIAL_DENSITY", nMaterialDensity, MaterialDensity);
  /* DESCRIPTION: Material thermal expansion coefficient */
  addDoubleListOption("MATERIAL_THERMAL_EXPANSION_COEFF", nMaterialThermalExpansion, MaterialThermalExpansion);
  /* DESCRIPTION: Temperature at which there is no stress from thermal expansion */
  addDoubleOption("MATERIAL_REFERENCE_TEMPERATURE", MaterialReferenceTemperature, 288.15);
  /* DESCRIPTION: Knowles B constant */
  addDoubleOption("KNOWLES_B", Knowles_B, 1.0);
  /* DESCRIPTION: Knowles N constant */
  addDoubleOption("KNOWLES_N", Knowles_N, 1.0);

  /*  DESCRIPTION: Include DE effects
  *  Options: NO, YES \ingroup Config */
  addBoolOption("DE_EFFECTS", DE_Effects, false);
  /*!\brief ELECTRIC_FIELD_CONST \n DESCRIPTION: Value of the Dielectric Elastomer constant */
  addDoubleListOption("ELECTRIC_FIELD_CONST", nElectric_Constant, Electric_Constant);
  /* DESCRIPTION: Modulus of the Electric Fields */
  addDoubleListOption("ELECTRIC_FIELD_MOD", nElectric_Field, Electric_Field_Mod);
  /* DESCRIPTION: Direction of the Electic Fields */
  addDoubleListOption("ELECTRIC_FIELD_DIR", nDim_Electric_Field, Electric_Field_Dir);

  /*!\brief DESIGN_VARIABLE_FEA
   *  \n DESCRIPTION: Design variable for FEA problems \n OPTIONS: See \link DVFEA_Map \endlink \n DEFAULT VENKATAKRISHNAN \ingroup Config */
  addEnumOption("DESIGN_VARIABLE_FEA", Kind_DV_FEA, DVFEA_Map, NODV_FEA);
  /*  DESCRIPTION: Consider a reference solution for the structure (optimization applications)
  *  Options: NO, YES \ingroup Config */
  addBoolOption("REFERENCE_GEOMETRY", RefGeom, false);
  /*!\brief REFERENCE_GEOMETRY_PENALTY\n DESCRIPTION: Penalty weight value for the objective function \ingroup Config*/
  addDoubleOption("REFERENCE_GEOMETRY_PENALTY", RefGeom_Penalty, 1E6);
  /*!\brief REFERENCE_GEOMETRY_FILENAME \n DESCRIPTION: Reference geometry filename \n Default: reference_geometry.dat \ingroup Config */
  addStringOption("REFERENCE_GEOMETRY_FILENAME", RefGeom_FEMFileName, string("reference_geometry.dat"));
  /*!\brief REFERENCE_GEOMETRY_FORMAT \n DESCRIPTION: Format of the reference geometry file \n OPTIONS: see \link Input_Ref_Map \endlink \n DEFAULT: SU2 \ingroup Config*/
  addEnumOption("REFERENCE_GEOMETRY_FORMAT", RefGeom_FileFormat, Input_Ref_Map, SU2_REF);
  /*!\brief REFERENCE_GEOMETRY_SURFACE\n DESCRIPTION: If true consider only the surfaces where loads are applied. \ingroup Config*/
  addBoolOption("REFERENCE_GEOMETRY_SURFACE", RefGeomSurf, false);

  /*!\brief REFERENCE_NODE\n  DESCRIPTION: Reference node for the structure (optimization applications) */
  addUnsignedLongOption("REFERENCE_NODE", refNodeID, 0);
  /*!\brief REFERENCE_NODE_DISPLACEMENT\n DESCRIPTION: Target displacement of the reference node \ingroup Config*/
  addDoubleListOption("REFERENCE_NODE_DISPLACEMENT", nDim_RefNode, RefNode_Displacement);
  /*!\brief REFERENCE_NODE_PENALTY\n DESCRIPTION: Penalty weight value for the objective function \ingroup Config*/
  addDoubleOption("REFERENCE_NODE_PENALTY", RefNode_Penalty, 1E3);

  /*!\brief TOTAL_DV_PENALTY\n DESCRIPTION: Penalty weight value to maintain the total sum of DV constant \ingroup Config*/
  addDoubleOption("TOTAL_DV_PENALTY", DV_Penalty, 0);

  /*!\brief STRESS_PENALTY_PARAM\n DESCRIPTION: Maximum allowed stress and KS exponent for structural optimization \ingroup Config*/
  addDoubleArrayOption("STRESS_PENALTY_PARAM", 2, StressPenaltyParam.data());

  /*!\brief REGIME_TYPE \n  DESCRIPTION: Geometric condition \n OPTIONS: see \link Struct_Map \endlink \ingroup Config*/
  addEnumOption("GEOMETRIC_CONDITIONS", Kind_Struct_Solver, Struct_Map, STRUCT_DEFORMATION::SMALL);
  /*!\brief REGIME_TYPE \n  DESCRIPTION: Material model \n OPTIONS: see \link Material_Map \endlink \ingroup Config*/
  addEnumOption("MATERIAL_MODEL", Kind_Material, Material_Map, STRUCT_MODEL::LINEAR_ELASTIC);
  /*!\brief REGIME_TYPE \n  DESCRIPTION: Compressibility of the material \n OPTIONS: see \link MatComp_Map \endlink \ingroup Config*/
  addEnumOption("MATERIAL_COMPRESSIBILITY", Kind_Material_Compress, MatComp_Map, STRUCT_COMPRESS::COMPRESSIBLE);

  /*  DESCRIPTION: Consider a prestretch in the structural domain
  *  Options: NO, YES \ingroup Config */
  addBoolOption("PRESTRETCH", Prestretch, false);
  /*!\brief PRESTRETCH_FILENAME \n DESCRIPTION: Filename to input for prestretching membranes \n Default: prestretch_file.dat \ingroup Config */
  addStringOption("PRESTRETCH_FILENAME", Prestretch_FEMFileName, string("prestretch_file.dat"));

  /* DESCRIPTION: Iterative method for non-linear structural analysis */
  addEnumOption("NONLINEAR_FEM_SOLUTION_METHOD", Kind_SpaceIteScheme_FEA, Space_Ite_Map_FEA, STRUCT_SPACE_ITE::NEWTON);
  /* DESCRIPTION: Formulation for bidimensional elasticity solver */
  addEnumOption("FORMULATION_ELASTICITY_2D", Kind_2DElasForm, ElasForm_2D, STRUCT_2DFORM::PLANE_STRAIN);
  /*  DESCRIPTION: Apply centrifugal forces
   *  Options: NO, YES \ingroup Config */
  addBoolOption("CENTRIFUGAL_FORCE", CentrifugalForce, false);
  /*  DESCRIPTION: Temporary: pseudo static analysis (no density in dynamic analysis)
   *  Options: NO, YES \ingroup Config */
  addBoolOption("PSEUDO_STATIC", PseudoStatic, false);
  /* DESCRIPTION: Parameter alpha for Newmark scheme (s) */
  addDoubleOption("NEWMARK_BETA", Newmark_beta, 0.25);
  /* DESCRIPTION: Parameter gamma for Newmark scheme (s) */
  addDoubleOption("NEWMARK_GAMMA", Newmark_gamma, 0.5);
  /* DESCRIPTION: Apply the load as a ramp */
  addBoolOption("RAMP_LOADING", Ramp_Load, false);
  /* DESCRIPTION: Time while the load is to be increased linearly */
  addDoubleOption("RAMP_TIME", Ramp_Time, 1.0);
  /* DESCRIPTION: Transfer method used for multiphysics problems */
  addEnumOption("DYNAMIC_LOAD_TRANSFER", Dynamic_LoadTransfer, Dyn_Transfer_Method_Map, POL_ORDER_1);

  /* DESCRIPTION: Newmark - Generalized alpha - coefficients */
  addDoubleListOption("TIME_INT_STRUCT_COEFFS", nIntCoeffs, Int_Coeffs);

  /*  DESCRIPTION: Apply dead loads. Options: NO, YES \ingroup Config */
  addBoolOption("INCREMENTAL_LOAD", IncrementalLoad, false);
  /* DESCRIPTION: Maximum number of increments of the  */
  addUnsignedLongOption("NUMBER_INCREMENTS", IncLoad_Nincrements, 10);

  inc_crit[0] = 0.0; inc_crit[1] = 0.0; inc_crit[2] = 0.0;
  /* DESCRIPTION: Definition of the  UTOL RTOL ETOL*/
  addDoubleArrayOption("INCREMENTAL_CRITERIA", 3, inc_crit);

  /* DESCRIPTION: Use of predictor */
  addBoolOption("PREDICTOR", Predictor, false);
  /* DESCRIPTION: Order of the predictor */
  addUnsignedShortOption("PREDICTOR_ORDER", Pred_Order, 0);

  /* DESCRIPTION: Topology optimization options */
  addBoolOption("TOPOLOGY_OPTIMIZATION", topology_optimization, false);
  addStringOption("TOPOL_OPTIM_OUTFILE", top_optim_output_file, string("element_derivatives.dat"));
  addDoubleOption("TOPOL_OPTIM_SIMP_EXPONENT", simp_exponent, 1.0);
  addDoubleOption("TOPOL_OPTIM_SIMP_MINSTIFF", simp_minimum_stiffness, 0.001);
  addEnumListOption("TOPOL_OPTIM_FILTER_KERNEL", top_optim_nKernel, top_optim_kernels, Filter_Kernel_Map);
  addDoubleListOption("TOPOL_OPTIM_FILTER_RADIUS", top_optim_nRadius, top_optim_filter_radius);
  addDoubleListOption("TOPOL_OPTIM_KERNEL_PARAM", top_optim_nKernelParams, top_optim_kernel_params);
  addUnsignedShortOption("TOPOL_OPTIM_SEARCH_LIMIT", top_optim_search_lim, 0);
  addEnumOption("TOPOL_OPTIM_PROJECTION_TYPE", top_optim_proj_type, Projection_Function_Map, ENUM_PROJECTION_FUNCTION::NONE);
  addDoubleOption("TOPOL_OPTIM_PROJECTION_PARAM", top_optim_proj_param, 0.0);

  /* CONFIG_CATEGORY: FSI solver */
  /*--- Options related to the FSI solver ---*/

  /* DESCRIPTION: ID of the region we want to compute the sensitivities using direct differentiation */
  addUnsignedShortOption("FEA_ID_DIRECTDIFF", nID_DV, 0);

  /* DESCRIPTION: Restart from a steady state (sets grid velocities to 0 when loading the restart). */
  addBoolOption("RESTART_STEADY_STATE", SteadyRestart, false);

  /*!\par CONFIG_CATEGORY: Multizone definition \ingroup Config*/
  /*--- Options related to multizone problems ---*/

  /* DESCRIPTION List of config files for each zone in a multizone setup with SOLVER=MULTIPHYSICS
   * Order here has to match the order in the meshfile if just one is used. */
  addStringListOption("CONFIG_LIST", nConfig_Files, Config_Filenames);

  /* DESCRIPTION: Determines if the multizone problem is solved for time-domain. */
  addBoolOption("TIME_DOMAIN", Time_Domain, false);
  /* DESCRIPTION: Number of outer iterations in the multizone problem. */
  addUnsignedLongOption("OUTER_ITER", nOuterIter, 1);
  /* DESCRIPTION: Number of inner iterations in each multizone block. */
  addUnsignedLongOption("INNER_ITER", nInnerIter, 1);
  /* DESCRIPTION: Number of time steps solved in the multizone problem. */
  addUnsignedLongOption("TIME_ITER", nTimeIter, 1);
  /* DESCRIPTION: Number of iterations in each single-zone block. */
  addUnsignedLongOption("ITER", nIter, 1000);
  /* DESCRIPTION: Restart iteration in the multizone problem. */
  addUnsignedLongOption("RESTART_ITER", Restart_Iter, 1);
  /* DESCRIPTION: Minimum error threshold for the linear solver for the implicit formulation */
  addDoubleOption("TIME_STEP", Time_Step, 0.0);
  /* DESCRIPTION: Total Physical Time for time-domain problems (s) */
  addDoubleOption("MAX_TIME", Max_Time, 1.0);
  /* DESCRIPTION: Determines if the special output is written out */
  addBoolOption("SPECIAL_OUTPUT", SpecialOutput, false);

  /* DESCRIPTION: Determines if the convergence history of each individual zone is written to screen */
  addBoolOption("WRT_ZONE_CONV", Wrt_ZoneConv, false);
  /* DESCRIPTION: Determines if the convergence history of each individual zone is written to file */
  addBoolOption("WRT_ZONE_HIST", Wrt_ZoneHist, false);

  /* DESCRIPTION: Determines if the forces breakdown is written out */
  addBoolOption("WRT_FORCES_BREAKDOWN", Wrt_ForcesBreakdown, false);


  /*!\par KIND_INTERPOLATION \n
   * DESCRIPTION: Type of interpolation to use for multi-zone problems. \n OPTIONS: see \link Interpolator_Map \endlink
   * Sets Kind_Interpolation \ingroup Config
   */
  addEnumOption("KIND_INTERPOLATION", Kind_Interpolation, Interpolator_Map, INTERFACE_INTERPOLATOR::NEAREST_NEIGHBOR);

  /*  DESCRIPTION: Use conservative approach for interpolating between meshes. */
  addBoolOption("CONSERVATIVE_INTERPOLATION", ConservativeInterpolation, true);

  addUnsignedShortOption("NUM_NEAREST_NEIGHBORS", NumNearestNeighbors, 1);

  /*!\par KIND_INTERPOLATION \n
   * DESCRIPTION: Type of radial basis function to use for radial basis function interpolation. \n OPTIONS: see \link RadialBasis_Map \endlink
   * Sets Kind_RadialBasis \ingroup Config
   */
  addEnumOption("KIND_RADIAL_BASIS_FUNCTION", Kind_RadialBasisFunction, RadialBasisFunction_Map, RADIAL_BASIS::WENDLAND_C2);

  /*  DESCRIPTION: Use polynomial term in radial basis function interpolation.
  *  Options: NO, YES \ingroup Config */
  addBoolOption("RADIAL_BASIS_FUNCTION_POLYNOMIAL_TERM", RadialBasisFunction_PolynomialOption, true);

  /* DESCRIPTION: Radius for radial basis function. */
  addDoubleOption("RADIAL_BASIS_FUNCTION_PARAMETER", RadialBasisFunction_Parameter, 1.0);

  /* DESCRIPTION: Tolerance to prune small coefficients from the RBF interpolation matrix. */
  addDoubleOption("RADIAL_BASIS_FUNCTION_PRUNE_TOLERANCE", RadialBasisFunction_PruneTol, 1e-6);

   /*!\par INLETINTERPOLATION \n
   * DESCRIPTION: Type of spanwise interpolation to use for the inlet face. \n OPTIONS: see \link Inlet_SpanwiseInterpolation_Map \endlink
   * Sets Kind_InletInterpolation \ingroup Config
   */
  addEnumOption("INLET_INTERPOLATION_FUNCTION",Kind_InletInterpolationFunction, Inlet_SpanwiseInterpolation_Map, INLET_SPANWISE_INTERP::NONE);

   /*!\par INLETINTERPOLATION \n
   * DESCRIPTION: Type of spanwise interpolation to use for the inlet face. \n OPTIONS: see \link Inlet_SpanwiseInterpolation_Map \endlink
   * Sets Kind_InletInterpolation \ingroup Config
   */
  addEnumOption("INLET_INTERPOLATION_DATA_TYPE", Kind_Inlet_InterpolationType, Inlet_SpanwiseInterpolationType_Map, INLET_INTERP_TYPE::VR_VTHETA);

  addBoolOption("PRINT_INLET_INTERPOLATED_DATA", PrintInlet_InterpolatedData, false);

  /* DESCRIPTION: Number of FSI iterations during which a ramp is applied */
  addUnsignedShortOption("RAMP_FSI_ITER", nIterFSI_Ramp, 2);
  /* DESCRIPTION: Aitken's static relaxation factor */
  addDoubleOption("STAT_RELAX_PARAMETER", AitkenStatRelax, 0.4);
  /* DESCRIPTION: Aitken's dynamic maximum relaxation factor for the first iteration */
  addDoubleOption("AITKEN_DYN_MAX_INITIAL", AitkenDynMaxInit, 0.5);
  /* DESCRIPTION: Aitken's dynamic minimum relaxation factor for the first iteration */
  addDoubleOption("AITKEN_DYN_MIN_INITIAL", AitkenDynMinInit, 0.5);
  /* DESCRIPTION: Kind of relaxation */
  addEnumOption("BGS_RELAXATION", Kind_BGS_RelaxMethod, AitkenForm_Map, BGS_RELAXATION::NONE);
  /* DESCRIPTION: Relaxation required */
  addBoolOption("RELAXATION", Relaxation, false);

  /*!\par CONFIG_CATEGORY: Radiation solver \ingroup Config*/
  /*--- Options related to the radiation solver ---*/

  /* DESCRIPTION: Type of radiation model */
  addEnumOption("RADIATION_MODEL", Kind_Radiation, Radiation_Map, RADIATION_MODEL::NONE);

  /* DESCRIPTION: Kind of initialization of the P1 model  */
  addEnumOption("P1_INITIALIZATION", Kind_P1_Init, P1_Init_Map, P1_INIT::TEMPERATURE);

  /* DESCRIPTION: Absorption coefficient */
  addDoubleOption("ABSORPTION_COEFF", Absorption_Coeff, 1.0);
  /* DESCRIPTION: Scattering coefficient */
  addDoubleOption("SCATTERING_COEFF", Scattering_Coeff, 0.0);

  /* DESCRIPTION: Apply a volumetric heat source as a source term (NO, YES) in the form of an ellipsoid*/
  addBoolOption("HEAT_SOURCE", HeatSource, false);
  /* DESCRIPTION: Value of the volumetric heat source */
  addDoubleOption("HEAT_SOURCE_VAL", ValHeatSource, 0.0);
  /* DESCRIPTION: Rotation of the volumetric heat source respect to Z axis */
  addDoubleOption("HEAT_SOURCE_ROTATION_Z", Heat_Source_Rot_Z, 0.0);
  /* DESCRIPTION: Position of heat source center (Heat_Source_Center_X, Heat_Source_Center_Y, Heat_Source_Center_Z) */
  hs_center[0] = 0.0; hs_center[1] = 0.0; hs_center[2] = 0.0;
  addDoubleArrayOption("HEAT_SOURCE_CENTER", 3, hs_center);
  /* DESCRIPTION: Vector of heat source radii (Heat_Source_Axes_A, Heat_Source_Axes_B, Heat_Source_Axes_C) */
  hs_axes[0] = 1.0; hs_axes[1] = 1.0; hs_axes[2] = 1.0;
  addDoubleArrayOption("HEAT_SOURCE_AXES", 3, hs_axes);

  /*!\brief MARKER_EMISSIVITY DESCRIPTION: Wall emissivity of the marker for radiation purposes \n
   * Format: ( marker, emissivity of the marker, ... ) \ingroup Config  */
  addStringDoubleListOption("MARKER_EMISSIVITY", nMarker_Emissivity, Marker_Emissivity, Wall_Emissivity);

  /* DESCRIPTION:  Courant-Friedrichs-Lewy condition of the finest grid in radiation solvers */
  addDoubleOption("CFL_NUMBER_RAD", CFL_Rad, 1.0);

  /*!\par CONFIG_CATEGORY: Heat solver \ingroup Config*/
  /*--- options related to the heat solver ---*/

  /* DESCRIPTION: CHT interface coupling methods */
  /*  Options: NO, YES \ingroup Config */
  addEnumOption("CHT_COUPLING_METHOD", Kind_CHT_Coupling, CHT_Coupling_Map, CHT_COUPLING::DIRECT_TEMPERATURE_ROBIN_HEATFLUX);

  /*!\par CONFIG_CATEGORY: Visualize Control Volumes \ingroup Config*/
  /*--- options related to visualizing control volumes ---*/

  /* DESCRIPTION: Node number for the CV to be visualized (tecplot) */
  addLongOption("VISUALIZE_CV", Visualize_CV, -1);

  /*!\par CONFIG_CATEGORY: Inverse design problem \ingroup Config*/
  /*--- options related to inverse design problem ---*/

  /* DESCRIPTION: Evaluate inverse design on the surface  */
  addBoolOption("INV_DESIGN_CP", InvDesign_Cp, false);

  /* DESCRIPTION: Evaluate inverse design on the surface  */
  addBoolOption("INV_DESIGN_HEATFLUX", InvDesign_HeatFlux, false);

  /*!\par CONFIG_CATEGORY: Unsupported options \ingroup Config*/
  /*--- Options that are experimental and not intended for general use ---*/

  /* DESCRIPTION: Write extra output (EXPERIMENTAL, NOT FOR GENERAL USE) */
  addBoolOption("EXTRA_OUTPUT", ExtraOutput, false);

  /* DESCRIPTION: Write extra heat output for a given heat solver zone */
  addLongOption("EXTRA_HEAT_ZONE_OUTPUT", ExtraHeatOutputZone, -1);

  /*--- options related to the FFD problem ---*/
  /*!\par CONFIG_CATEGORY:FFD point inversion \ingroup Config*/

  /* DESCRIPTION: Fix I plane */
  addShortListOption("FFD_FIX_I", nFFD_Fix_IDir, FFD_Fix_IDir);

  /* DESCRIPTION: Fix J plane */
  addShortListOption("FFD_FIX_J", nFFD_Fix_JDir, FFD_Fix_JDir);

  /* DESCRIPTION: Fix K plane */
  addShortListOption("FFD_FIX_K", nFFD_Fix_KDir, FFD_Fix_KDir);

  /* DESCRIPTION: FFD symmetry plane (j=0) */
  addBoolOption("FFD_SYMMETRY_PLANE", FFD_Symmetry_Plane, false);

  /* DESCRIPTION: Define different coordinates systems for the FFD */
  addEnumOption("FFD_COORD_SYSTEM", FFD_CoordSystem, CoordSystem_Map, CARTESIAN);

  /* DESCRIPTION: Axis information for the spherical and cylindrical coord system */
  ffd_axis[0] = 0.0; ffd_axis[1] = 0.0; ffd_axis[2] =0.0;
  addDoubleArrayOption("FFD_AXIS", 3, ffd_axis);

  /* DESCRIPTION: Number of total iterations in the FFD point inversion */
  addUnsignedShortOption("FFD_ITERATIONS", nFFD_Iter, 500);

  /* DESCRIPTION: Free surface damping coefficient */
  addDoubleOption("FFD_TOLERANCE", FFD_Tol, 1E-10);

  /* DESCRIPTION: Procedure to prevent self-intersections within the FFD box based on Jacobian determinant */
  addBoolOption("FFD_INTPREV", FFD_IntPrev, NO);

  /* DESCRIPTION: Number of total iterations in the convexity check procedure */
  addUnsignedShortOption("FFD_INTPREV_ITER", FFD_IntPrev_MaxIter, 10);

  /* DESCRIPTION: Recursion depth in the FFD self-intersection prevention */
  addUnsignedShortOption("FFD_INTPREV_DEPTH", FFD_IntPrev_MaxDepth, 3);

  /* DESCRIPTION: Convexity check on all mesh elements */
  addBoolOption("CONVEXITY_CHECK", ConvexityCheck, NO);

  /* DESCRIPTION: Number of total iterations in the convexity check procedure */
  addUnsignedShortOption("CONVEXITY_CHECK_ITER", ConvexityCheck_MaxIter, 10);

  /* DESCRIPTION: Recursion depth in the FFD self-intersection prevention */
  addUnsignedShortOption("CONVEXITY_CHECK_DEPTH", ConvexityCheck_MaxDepth, 3);

  /* DESCRIPTION: Definition of the FFD boxes */
  addFFDDefOption("FFD_DEFINITION", nFFDBox, CoordFFDBox, TagFFDBox);

  /* DESCRIPTION: Definition of the FFD boxes */
  addFFDDegreeOption("FFD_DEGREE", nFFDBox, DegreeFFDBox);

  /* DESCRIPTION: Surface continuity at the intersection with the FFD */
  addEnumOption("FFD_CONTINUITY", FFD_Continuity, Continuity_Map, DERIVATIVE_2ND);

  /* DESCRIPTION: Kind of blending for the FFD definition */
  addEnumOption("FFD_BLENDING", FFD_Blending, Blending_Map, BEZIER );

  /* DESCRIPTION: Order of the BSplines for BSpline Blending function */
  ffd_coeff[0] = 2; ffd_coeff[1] = 2; ffd_coeff[2] = 2;
  addDoubleArrayOption("FFD_BSPLINE_ORDER", 3, ffd_coeff);

  /*--- Options for the automatic differentiation methods ---*/
  /*!\par CONFIG_CATEGORY: Automatic Differentation options\ingroup Config*/

  /* DESCRIPTION: Direct differentiation mode (forward) */
  addEnumOption("DIRECT_DIFF", DirectDiff, DirectDiff_Var_Map, NO_DERIVATIVE);

  /* DESCRIPTION: Automatic differentiation mode (reverse) */
  addBoolOption("AUTO_DIFF", AD_Mode, NO);

  /* DESCRIPTION: Preaccumulation in the AD mode. */
  addBoolOption("PREACC", AD_Preaccumulation, YES);

  /*--- options that are used in the python optimization scripts. These have no effect on the c++ toolsuite ---*/
  /*!\par CONFIG_CATEGORY:Python Options\ingroup Config*/

  /* DESCRIPTION: Gradient method */
  addPythonOption("GRADIENT_METHOD");

  /* DESCRIPTION: Geometrical Parameter */
  addPythonOption("GEO_PARAM");

  /* DESCRIPTION: Setup for design variables */
  addPythonOption("DEFINITION_DV");

  /* DESCRIPTION: Maximum number of iterations */
  addPythonOption("OPT_ITERATIONS");

  /* DESCRIPTION: Requested accuracy */
  addPythonOption("OPT_ACCURACY");

  /*!\brief OPT_COMBINE_OBJECTIVE
   *  \n DESCRIPTION: Flag specifying whether to internally combine a multi-objective function or treat separately */
  addPythonOption("OPT_COMBINE_OBJECTIVE");

  /* DESCRIPTION: Current value of the design variables */
  addPythonOption("DV_VALUE_NEW");

  /* DESCRIPTION: Previous value of the design variables */
  addPythonOption("DV_VALUE_OLD");

  /* DESCRIPTION: Number of partitions of the mesh */
  addPythonOption("NUMBER_PART");

  /* DESCRIPTION: Optimization objective function with optional scaling factor*/
  addPythonOption("OPT_OBJECTIVE");

  /* DESCRIPTION: Optimization constraint functions with optional scaling factor */
  addPythonOption("OPT_CONSTRAINT");

  /* DESCRIPTION: Finite different step for gradient estimation */
  addPythonOption("FIN_DIFF_STEP");

  /* DESCRIPTION: Verbosity of the python scripts to Stdout */
  addPythonOption("CONSOLE");

  /* DESCRIPTION: Flag specifying if the mesh was decomposed */
  addPythonOption("DECOMPOSED");

  /* DESCRIPTION: Optimization gradient factor */
  addPythonOption("OPT_GRADIENT_FACTOR");

  /* DESCRIPTION: Upper bound for the optimizer */
  addPythonOption("OPT_BOUND_UPPER");

  /* DESCRIPTION: Lower bound for the optimizer */
  addPythonOption("OPT_BOUND_LOWER");

  /* DESCRIPTION: Number of zones of the problem */
  addPythonOption("NZONES");

  /* DESCRIPTION: ParMETIS load balancing tolerance */
  addDoubleOption("PARMETIS_TOLERANCE", ParMETIS_tolerance, 0.02);

  /* DESCRIPTION: ParMETIS load balancing weight for points */
  addLongOption("PARMETIS_POINT_WEIGHT", ParMETIS_pointWgt, 0);

  /* DESCRIPTION: ParMETIS load balancing weight for edges (equiv. to neighbors) */
  addLongOption("PARMETIS_EDGE_WEIGHT", ParMETIS_edgeWgt, 1);

  /*--- options that are used in the Hybrid RANS/LES Simulations  ---*/
  /*!\par CONFIG_CATEGORY:Hybrid_RANSLES Options\ingroup Config*/

  /* DESCRIPTION: Starting Iteration for windowing approach */
  addUnsignedLongOption("WINDOW_START_ITER", StartWindowIteration, 0);

  /* DESCRIPTION: Window (weight) function for the cost-functional in the reverse sweep */
  addEnumOption("WINDOW_FUNCTION", Kind_WindowFct, Window_Map, WINDOW_FUNCTION::SQUARE);

  /* DESCRIPTION: DES Constant */
  addDoubleOption("DES_CONST", Const_DES, 0.65);

  /* DESCRIPTION: Specify Hybrid RANS/LES model */
  addEnumOption("HYBRID_RANSLES", Kind_HybridRANSLES, HybridRANSLES_Map, NO_HYBRIDRANSLES);

  /* DESCRIPTION: Roe with low dissipation for unsteady flows */
  addEnumOption("ROE_LOW_DISSIPATION", Kind_RoeLowDiss, RoeLowDiss_Map, NO_ROELOWDISS);

  /* DESCRIPTION: Compute Average for unsteady simulations */
  addBoolOption("COMPUTE_AVERAGE", Compute_Average, false);

  /* DESCRIPTION: Multipoint design Mach number*/
  addPythonOption("MULTIPOINT_MACH_NUMBER");

  /* DESCRIPTION: Multipoint design Weight */
  addPythonOption("MULTIPOINT_WEIGHT");

  /* DESCRIPTION: Multipoint design Angle of Attack */
  addPythonOption("MULTIPOINT_AOA");

  /* DESCRIPTION: Multipoint design Sideslip angle */
  addPythonOption("MULTIPOINT_SIDESLIP_ANGLE");

  /* DESCRIPTION: Multipoint design target CL*/
  addPythonOption("MULTIPOINT_TARGET_CL");

  /* DESCRIPTION: Multipoint design Reynolds number */
  addPythonOption("MULTIPOINT_REYNOLDS_NUMBER");

  /* DESCRIPTION: Multipoint design freestream temperature */
  addPythonOption("MULTIPOINT_FREESTREAM_TEMPERATURE");

  /* DESCRIPTION: Multipoint design freestream pressure */
  addPythonOption("MULTIPOINT_FREESTREAM_PRESSURE");

  /* DESCRIPTION: Multipoint design for outlet quantities (varying back pressure or mass flow operating points). */
  addPythonOption("MULTIPOINT_OUTLET_VALUE");

  /* DESCRIPTION: Multipoint mesh filenames, if using different meshes for each point */
  addPythonOption("MULTIPOINT_MESH_FILENAME");

  /*--- options that are used for the output ---*/
  /*!\par CONFIG_CATEGORY:Output Options\ingroup Config*/

  /* DESCRIPTION: Type of screen output */
  addStringListOption("SCREEN_OUTPUT", nScreenOutput, ScreenOutput);
  /* DESCRIPTION: Type of output printed to the history file */
  addStringListOption("HISTORY_OUTPUT", nHistoryOutput, HistoryOutput);
  /* DESCRIPTION: Type of output printed to the volume solution file */
  addStringListOption("VOLUME_OUTPUT", nVolumeOutput, VolumeOutput);

  /* DESCRIPTION: History writing frequency (INNER_ITER) */
  addUnsignedLongOption("HISTORY_WRT_FREQ_INNER", HistoryWrtFreq[2], 1);
  /* DESCRIPTION: History writing frequency (OUTER_ITER) */
  addUnsignedLongOption("HISTORY_WRT_FREQ_OUTER", HistoryWrtFreq[1], 1);
  /* DESCRIPTION: History writing frequency (TIME_ITER) */
  addUnsignedLongOption("HISTORY_WRT_FREQ_TIME", HistoryWrtFreq[0], 1);

  /* DESCRIPTION: Screen writing frequency (INNER_ITER) */
  addUnsignedLongOption("SCREEN_WRT_FREQ_INNER", ScreenWrtFreq[2], 1);
  /* DESCRIPTION: Screen writing frequency (OUTER_ITER) */
  addUnsignedLongOption("SCREEN_WRT_FREQ_OUTER", ScreenWrtFreq[1], 1);
  /* DESCRIPTION: Screen writing frequency (TIME_ITER) */
  addUnsignedLongOption("SCREEN_WRT_FREQ_TIME", ScreenWrtFreq[0], 1);
  /* DESCRIPTION: list of writing frequencies for each file type (length same as nVolumeOutputFiles) */
  addULongListOption("OUTPUT_WRT_FREQ", nVolumeOutputFrequencies, VolumeOutputFrequencies);

  /* DESCRIPTION: Volume solution files */
  addEnumListOption("OUTPUT_FILES", nVolumeOutputFiles, VolumeOutputFiles, Output_Map);

  /* DESCRIPTION: Parameter to perturb eigenvalues */
  addDoubleOption("UQ_DELTA_B", uq_delta_b, 1.0);

  /* DESCRIPTION: Parameter to determine kind of perturbation */
  addUnsignedShortOption("UQ_COMPONENT", eig_val_comp, 1);

  /* DESCRIPTION: Parameter to perturb eigenvalues */
  addDoubleOption("UQ_URLX", uq_urlx, 0.1);

  /* DESCRIPTION: Permuting eigenvectors for UQ analysis */
  addBoolOption("UQ_PERMUTE", uq_permute, false);

  /* DESCRIPTION: Number of calls to 'Build' that trigger re-factorization (0 means only once). */
  addUnsignedLongOption("PASTIX_FACTORIZATION_FREQUENCY", pastix_fact_freq, 1);

  /* DESCRIPTION: 0 - Quiet, 1 - During factorization and cleanup, 2 - Even more detail. */
  addUnsignedShortOption("PASTIX_VERBOSITY_LEVEL", pastix_verb_lvl, 0);

  /* DESCRIPTION: Level of fill for PaStiX incomplete LU factorization. */
  addUnsignedShortOption("PASTIX_FILL_LEVEL", pastix_fill_lvl, 1);

  /* DESCRIPTION: Size of the edge groups colored for thread parallel edge loops (0 forces the reducer strategy). */
  addUnsignedLongOption("EDGE_COLORING_GROUP_SIZE", edgeColorGroupSize, 512);

  /* DESCRIPTION: Allow fallback to smaller edge color group sizes for the discrete adjoint and allow more colors. */
  addBoolOption("EDGE_COLORING_RELAX_DISC_ADJ", edgeColoringRelaxDiscAdj, true);

  /*--- options that are used for libROM ---*/
  /*!\par CONFIG_CATEGORY:libROM options \ingroup Config*/

  /*!\brief SAVE_LIBROM \n DESCRIPTION: Flag for saving data with libROM. */
  addBoolOption("SAVE_LIBROM", libROM, false);

  /*!\brief LIBROM_BASE_FILENAME \n DESCRIPTION: Output base file name for libROM   \ingroup Config*/
  addStringOption("LIBROM_BASE_FILENAME", libROMbase_FileName, string("su2"));

  /*!\brief BASIS_GENERATION \n DESCRIPTION: Flag for saving data with libROM. */
  addEnumOption("BASIS_GENERATION", POD_Basis_Gen, POD_Map, POD_KIND::STATIC);

  /*!\brief MAX_BASIS_DIM \n DESCRIPTION: Maximum number of basis vectors.*/
  addUnsignedShortOption("MAX_BASIS_DIM", maxBasisDim, 100);

  /*!\brief ROM_SAVE_FREQ \n DESCRIPTION: How often to save snapshots for unsteady problems.*/
  addUnsignedShortOption("ROM_SAVE_FREQ", rom_save_freq, 1);

  /* END_CONFIG_OPTIONS */

}

void CConfig::SetConfig_Parsing(char case_filename[MAX_STRING_SIZE]) {

  ifstream case_file;

  /*--- Read the configuration file ---*/

  case_file.open(case_filename, ios::in);

  if (case_file.fail()) {
    SU2_MPI::Error("The configuration file (.cfg) is missing!!", CURRENT_FUNCTION);
  }

  SetConfig_Parsing(case_file);

  case_file.close();

}

void CConfig::SetConfig_Parsing(istream& config_buffer){

  string text_line, option_name;
  vector<string> option_value;

  string errorString;

  const int max_err_count = 30; // Maximum number of errors to print before stopping
  int err_count = 0;  // How many errors have we found in the config file
  int line_count = 1;

  map<string, bool> included_options;

  /*--- Parse the configuration file and set the options ---*/

  while (getline (config_buffer, text_line)) {

    if (err_count >= max_err_count) {
      errorString.append("Too many errors, stopping parse.");
      break;
    }

     PrintingToolbox::trim(text_line);

    /*--- Check if there is a line continuation character at the
     * end of the current line or somewhere in between (the rest is ignored then).
     * If yes, read until there is a line without one or an empty line.
     * If there is a statement after a cont. char
     * throw an error. ---*/

     if (!text_line.empty() && (text_line.front() != '%')){
       while (text_line.back() == '\\' ||
              (PrintingToolbox::split(text_line, '\\').size() > 1)){
         string tmp;
         getline (config_buffer, tmp);
         line_count++;
         if (tmp.find_first_of('=') != string::npos){
           errorString.append("Line " + to_string(line_count)  + ": Statement found after continuation character.\n");
         }
         PrintingToolbox::trim(tmp);
         if (tmp.front() != '%'){
           text_line = PrintingToolbox::split(text_line, '\\')[0];
           text_line += " " + tmp;
         }
       }
     }

    if (TokenizeString(text_line, option_name, option_value)) {
      /*--- See if it's a python option ---*/

      if (option_map.find(option_name) == option_map.end()) {
          string newString;
          newString.append("Line " + to_string(line_count)  + " " + option_name);
          newString.append(": invalid option name");
          newString.append(". Check current SU2 options in config_template.cfg.");
          newString.append("\n");
          if (!option_name.compare("SINGLEZONE_DRIVER"))
            newString.append("Option SINGLEZONE_DRIVER is deprecated, it does not have a replacement.\n\n");
          else if (!option_name.compare("DYN_TIMESTEP"))
            newString.append("DYN_TIMESTEP is deprecated. Use TIME_STEP instead.\n\n");
          else if (!option_name.compare("DYN_TIME"))
            newString.append("DYN_TIME is deprecated. Use MAX_TIME instead.\n\n");
          else if (!option_name.compare("DYNAMIC_ANALYSIS"))
            newString.append("DYNAMIC_ANALYSIS is deprecated. Use TIME_DOMAIN instead.\n\n");
          else if (!option_name.compare("SPECIES_USE_STRONG_BC"))
            newString.append("SPECIES_USE_STRONG_BC is deprecated. Use MARKER_SPECIES_STRONG_BC= (marker1, ...) instead.\n\n");
          else if (!option_name.compare("DEAD_LOAD"))
            newString.append("DEAD_LOAD is deprecated. Use GRAVITY_FORCE or BODY_FORCE instead.\n\n");
          else {
            /*--- Find the most likely candidate for the unrecognized option, based on the length
             of start and end character sequences shared by candidates and the option. ---*/
            auto countMatchChars = [&option_name](const string& candidate) {
              const size_t sz1 = option_name.size(), sz2 = candidate.size();
              size_t nMatch = 0;
              for (size_t i=0; i<min(sz1,sz2); ++i) {
                if (option_name[i] == candidate[i]) nMatch++;
                else break;
              }
              for (size_t i=0; i<min(sz1,sz2); ++i) {
                if (option_name[sz1-1-i] == candidate[sz2-1-i]) nMatch++;
                else break;
              }
              return nMatch;
            };
            string match;
            size_t maxScore = 0;
            for (auto& candidate : option_map) {
              auto score = countMatchChars(candidate.first);
              if (score > maxScore) {
                maxScore = score;
                match = candidate.first;
              }
            }
            newString.append("Did you mean ");
            newString.append(match);
            newString.append("?\n");
          }
          errorString.append(newString);
          err_count++;
          line_count++;
        continue;
      }

      /*--- Option exists, check if the option has already been in the config file ---*/

      if (included_options.find(option_name) != included_options.end()) {
        string newString;
        newString.append("Line " + to_string(line_count)  + " " + option_name);
        newString.append(": option appears twice");
        newString.append("\n");
        errorString.append(newString);
        err_count++;
        line_count++;
        continue;
      }

      /*--- New found option. Add it to the map, and delete from all options ---*/

      included_options.insert(pair<string, bool>(option_name, true));
      all_options.erase(option_name);

      /*--- Set the value and check error ---*/

      string out = option_map[option_name]->SetValue(option_value);
      if (out.compare("") != 0) {
        /*--- valid option, but deprecated value ---*/
        if (!option_name.compare("KIND_TURB_MODEL")) {
          if (option_value[0] == "SST_SUST")
            errorString.append("Option KIND_TURB_MODEL=SST_SUST is deprecated. Use KIND_TURB_MODEL=SST, SST_OPTIONS=SUSTAINING instead.\n");
          else if (option_value[0] == "SA_NEG")
            errorString.append("Option KIND_TURB_MODEL=SA_NEG is deprecated. Use KIND_TURB_MODEL=SA, SA_OPTIONS=NEGATIVE instead.\n");
          else if (option_value[0] == "SA_E")
            errorString.append("Option KIND_TURB_MODEL=SA_E is deprecated. Use KIND_TURB_MODEL=SA, SA_OPTIONS=EDWARDS instead.\n");
          else if (option_value[0] == "SA_COMP")
            errorString.append("Option KIND_TURB_MODEL=SA_COMP is deprecated. Use KIND_TURB_MODEL=SA, SA_OPTIONS=COMPRESSIBILITY instead.\n");
          else if (option_value[0] == "SA_E_COMP")
            errorString.append("Option KIND_TURB_MODEL=SA_E_COMP is deprecated. Use KIND_TURB_MODEL=SA, SA_OPTIONS=EDWARDS,COMPRESSIBILITY instead.\n");
        } else if (!option_name.compare("KIND_TRANS_MODEL")) {
          if (option_value[0] == "BC")
            errorString.append("Option KIND_TRANS_MODEL=BC is deprecated. Use KIND_TURB_MODEL=SA, SA_OPTIONS=BCM instead.\n");
        }
        errorString.append(out);
        errorString.append("\n");
        err_count++;
      }
    }
    line_count++;
  }

  /*--- See if there were any errors parsing the config file ---*/

  if (!errorString.empty()) {
    SU2_MPI::Error(errorString, CURRENT_FUNCTION);
  }
}

void CConfig::SetDefaultFromConfig(CConfig *config){

  map<string, bool> noInheritance = {{"SCREEN_OUTPUT", true},{"HISTORY_OUTPUT", true}};

  map<string, bool>::iterator iter = all_options.begin(), curr_iter;

  while (iter != all_options.end()){
    curr_iter = iter++;
    if (!config->option_map[curr_iter->first]->GetValue().empty() && !noInheritance[curr_iter->first]){
      option_map[curr_iter->first]->SetValue(config->option_map[curr_iter->first]->GetValue());
      all_options.erase(curr_iter);
    }
  }
}

void CConfig::SetDefault(){

  /*--- Set the default values for all of the options that weren't set ---*/

  for (auto iter = all_options.begin(); iter != all_options.end(); ++iter) {
    if (option_map[iter->first]->GetValue().empty())
      option_map[iter->first]->SetDefault();
  }
}

bool CConfig::SetRunTime_Parsing(char case_filename[MAX_STRING_SIZE]) {
  string text_line, option_name;
  ifstream case_file;
  vector<string> option_value;

  /*--- Read the configuration file ---*/

  case_file.open(case_filename, ios::in);

  if (case_file.fail()) { return false; }

  string errorString;

  int err_count = 0;  // How many errors have we found in the config file
  const int max_err_count = 30; // Maximum number of errors to print before stopping

  map<string, bool> included_options;

  /*--- Parse the configuration file and set the options ---*/

  while (getline (case_file, text_line)) {

    if (err_count >= max_err_count) {
      errorString.append("Too many errors, stopping parse.");
      break;
    }

    if (TokenizeString(text_line, option_name, option_value)) {

      if (option_map.find(option_name) == option_map.end()) {

        /*--- See if it's a python option ---*/

        string newString;
        newString.append(option_name);
        newString.append(": invalid option name");
        newString.append("\n");
        errorString.append(newString);
        err_count++;
        continue;
      }

      /*--- Option exists, check if the option has already been in the config file ---*/

      if (included_options.find(option_name) != included_options.end()) {
        string newString;
        newString.append(option_name);
        newString.append(": option appears twice");
        newString.append("\n");
        errorString.append(newString);
        err_count++;
        continue;
      }

      /*--- New found option. Add it to the map, and delete from all options ---*/

      included_options.insert(pair<string, bool>(option_name, true));
      all_options.erase(option_name);

      /*--- Set the value and check error ---*/

      string out = option_map[option_name]->SetValue(option_value);
      if (out.compare("") != 0) {
        errorString.append(out);
        errorString.append("\n");
        err_count++;
      }

    }
  }

  /*--- Set the default values for all of the options that weren't set ---*/

  for (auto iter = all_options.begin(); iter != all_options.end(); ++iter) {
    option_map[iter->first]->SetDefault();
  }

  /*--- See if there were any errors parsing the runtime file ---*/

  if (!errorString.empty()) {
    SU2_MPI::Error(errorString, CURRENT_FUNCTION);
  }

  case_file.close();

  return true;

}

void CConfig::SetHeader(SU2_COMPONENT val_software) const{

  if ((iZone == 0) && (rank == MASTER_NODE)) {
    cout << "\n";
    cout << "-------------------------------------------------------------------------\n";
    cout << "|    ___ _   _ ___                                                      |\n";
    cout << "|   / __| | | |_  )   Release 8.1.0 \"Harrier\"                           |\n";
    cout << "|   \\__ \\ |_| |/ /                                                      |\n";
    switch (val_software) {
    case SU2_COMPONENT::SU2_CFD: cout << "|   |___/\\___//___|   Suite (Computational Fluid Dynamics Code)         |\n"; break;
    case SU2_COMPONENT::SU2_DEF: cout << "|   |___/\\___//___|   Suite (Mesh Deformation Code)                     |\n"; break;
    case SU2_COMPONENT::SU2_DOT: cout << "|   |___/\\___//___|   Suite (Gradient Projection Code)                  |\n"; break;
    case SU2_COMPONENT::SU2_GEO: cout << "|   |___/\\___//___|   Suite (Geometry Definition Code)                  |\n"; break;
    case SU2_COMPONENT::SU2_SOL: cout << "|   |___/\\___//___|   Suite (Solution Exporting Code)                   |\n"; break;
    }
    cout << "|                                                                       |\n";
    cout << "-------------------------------------------------------------------------\n";
    cout << "| SU2 Project Website: https://su2code.github.io                        |\n";
    cout << "|                                                                       |\n";
    cout << "| The SU2 Project is maintained by the SU2 Foundation                   |\n";
    cout << "| (http://su2foundation.org)                                            |\n";
    cout << "-------------------------------------------------------------------------\n";
    cout << "| Copyright 2012-2024, SU2 Contributors                                 |\n";
    cout << "|                                                                       |\n";
    cout << "| SU2 is free software; you can redistribute it and/or                  |\n";
    cout << "| modify it under the terms of the GNU Lesser General Public            |\n";
    cout << "| License as published by the Free Software Foundation; either          |\n";
    cout << "| version 2.1 of the License, or (at your option) any later version.    |\n";
    cout << "|                                                                       |\n";
    cout << "| SU2 is distributed in the hope that it will be useful,                |\n";
    cout << "| but WITHOUT ANY WARRANTY; without even the implied warranty of        |\n";
    cout << "| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |\n";
    cout << "| Lesser General Public License for more details.                       |\n";
    cout << "|                                                                       |\n";
    cout << "| You should have received a copy of the GNU Lesser General Public      |\n";
    cout << "| License along with SU2. If not, see <http://www.gnu.org/licenses/>.   |\n";
    cout << "-------------------------------------------------------------------------" << endl;
  }

}

void CConfig::SetnZone(){

  /*--- Just as a clarification --- */

  if (Multizone_Problem == NO && Kind_Solver != MAIN_SOLVER::MULTIPHYSICS){
    nZone = 1;
  }

  if (Kind_Solver == MAIN_SOLVER::MULTIPHYSICS){
    Multizone_Problem = YES;
    if (nConfig_Files == 0){
      SU2_MPI::Error("CONFIG_LIST must be provided if PHYSICAL_PROBLEM=MULTIPHYSICS", CURRENT_FUNCTION);
    }
  }

  if (Multizone_Problem == YES){

    /*--- Some basic multizone checks ---*/

    if (nMarker_ZoneInterface % 2 != 0){
      SU2_MPI::Error("Number of markers in MARKER_ZONE_INTERFACE must be a multiple of 2", CURRENT_FUNCTION);
    }

    if (Multizone_Mesh){

      /*--- Get the number of zones from the mesh file --- */

      nZone = GetnZone(Mesh_FileName, Mesh_FileFormat);

      /*--- If config list is set, make sure number matches number of zones in mesh file --- */

      if (nConfig_Files != 0 && (nZone != nConfig_Files)){
        SU2_MPI::Error("Number of CONFIG_LIST must match number of zones in mesh file.", CURRENT_FUNCTION);
      }
    } else {

      /*--- Number of zones is determined from the number of config files provided --- */

      if (nConfig_Files == 0){
        SU2_MPI::Error("If MULTIZONE_MESH is set to YES, you must provide a list of config files using CONFIG_LIST option", CURRENT_FUNCTION);
      }
      nZone = nConfig_Files;

    }

    /*--- Check if subconfig files exist --- */

    if (nConfig_Files != 0){
      for (unsigned short iConfig = 0; iConfig < nConfig_Files; iConfig++){
        ifstream f(Config_Filenames[iConfig].c_str());
        if (!f.good()){
          SU2_MPI::Error("Config file " + Config_Filenames[iConfig] + " defined in CONFIG_FILES does not exist", CURRENT_FUNCTION);
        }
      }
    }

  }

}

void CConfig::SetPostprocessing(SU2_COMPONENT val_software, unsigned short val_izone, unsigned short val_nDim) {

  unsigned short iCFL, iMarker;
  bool ideal_gas = ((Kind_FluidModel == STANDARD_AIR) ||
                    (Kind_FluidModel == IDEAL_GAS) ||
                    (Kind_FluidModel == INC_IDEAL_GAS) ||
                    (Kind_FluidModel == FLUID_MIXTURE) ||
                    (Kind_FluidModel == FLUID_FLAMELET) ||
                    (Kind_FluidModel == INC_IDEAL_GAS_POLY) ||
                    (Kind_FluidModel == CONSTANT_DENSITY));
  bool noneq_gas = ((Kind_FluidModel == MUTATIONPP) ||
                    (Kind_FluidModel == SU2_NONEQ));
  bool standard_air = ((Kind_FluidModel == STANDARD_AIR));
  bool nemo = GetNEMOProblem();

  if (nZone > 1){
    Multizone_Problem = YES;
  }

  /*--- Set the default output files ---*/
  if (!OptionIsSet("OUTPUT_FILES")){
    nVolumeOutputFiles = 3;
    VolumeOutputFiles = new OUTPUT_TYPE[nVolumeOutputFiles];
    VolumeOutputFiles[0] = OUTPUT_TYPE::RESTART_BINARY;
    VolumeOutputFiles[1] = OUTPUT_TYPE::PARAVIEW_XML;
    VolumeOutputFiles[2] = OUTPUT_TYPE::SURFACE_PARAVIEW_XML;
  }

  /*--- Set the default output frequencies ---*/
  if (!OptionIsSet("OUTPUT_WRT_FREQ")){
    nVolumeOutputFrequencies = nVolumeOutputFiles;
    VolumeOutputFrequencies = new unsigned long [nVolumeOutputFrequencies];

    /*---  Using default frequency of 250 for all files when steady, and 1 for unsteady. ---*/
    for (auto iVolumeFreq = 0; iVolumeFreq < nVolumeOutputFrequencies; iVolumeFreq++){
      VolumeOutputFrequencies[iVolumeFreq] = Time_Domain ? 1 : 250;
    }
  } else if (nVolumeOutputFrequencies < nVolumeOutputFiles) {
    /*--- If there are fewer frequencies than files, repeat the last frequency.
     *    This is useful to define 1 frequency for the restart file and 1 frequency for all the visualization files.  ---*/
    auto* newFrequencies = new unsigned long[nVolumeOutputFiles];
    for (unsigned short i = 0; i < nVolumeOutputFrequencies; ++i) {
      newFrequencies[i] = VolumeOutputFrequencies[i];
    }
    for (auto i = nVolumeOutputFrequencies; i < nVolumeOutputFiles; ++i) {
      newFrequencies[i] = newFrequencies[i-1];
    }
    delete [] VolumeOutputFrequencies;
    VolumeOutputFrequencies = newFrequencies;
    nVolumeOutputFrequencies = nVolumeOutputFiles;
  }

  /*--- Check if SU2 was build with TecIO support, as that is required for Tecplot Binary output. ---*/
#ifndef HAVE_TECIO
  for (unsigned short iVolumeFile = 0; iVolumeFile < nVolumeOutputFiles; iVolumeFile++){
    if (VolumeOutputFiles[iVolumeFile] == OUTPUT_TYPE::TECPLOT_BINARY ||
        VolumeOutputFiles[iVolumeFile] == OUTPUT_TYPE::SURFACE_TECPLOT_BINARY) {
      SU2_MPI::Error(string("Tecplot binary file requested in option OUTPUT_FILES but SU2 was built without TecIO support.\n"), CURRENT_FUNCTION);
    }
  }
#endif

  /*--- Check if SU2 was build with CGNS support, as that is required for CGNS output. ---*/
#ifndef HAVE_CGNS
  for (unsigned short iVolumeFile = 0; iVolumeFile < nVolumeOutputFiles; iVolumeFile++) {
    if (VolumeOutputFiles[iVolumeFile] == OUTPUT_TYPE::CGNS ||
        VolumeOutputFiles[iVolumeFile] == OUTPUT_TYPE::SURFACE_CGNS) {
      SU2_MPI::Error(string("CGNS file requested in option OUTPUT_FILES but SU2 was built without CGNS support.\n"),CURRENT_FUNCTION);
    }
  }
#endif

  /*--- Check if CoolProp is used with non-dimensionalization. ---*/
  if (Kind_FluidModel == COOLPROP && Ref_NonDim != DIMENSIONAL) {
    SU2_MPI::Error("CoolProp can not be used with non-dimensionalization.", CURRENT_FUNCTION);
  }

  /*--- STL_BINARY output not implemented yet, but already a value in option_structure.hpp---*/
  for (unsigned short iVolumeFile = 0; iVolumeFile < nVolumeOutputFiles; iVolumeFile++) {
    if (VolumeOutputFiles[iVolumeFile] == OUTPUT_TYPE::STL_BINARY){
      SU2_MPI::Error(string("OUTPUT_FILES: 'STL_BINARY' output not implemented. Use 'STL' for ASCII output.\n"), CURRENT_FUNCTION);
    }
    if (val_nDim == 2 && (VolumeOutputFiles[iVolumeFile] == OUTPUT_TYPE::STL_ASCII || VolumeOutputFiles[iVolumeFile] == OUTPUT_TYPE::STL_BINARY)) {
      SU2_MPI::Error(string("OUTPUT_FILES: 'STL(_BINARY)' output only reasonable for 3D cases.\n"), CURRENT_FUNCTION);
    }
  }

  /*--- Check if MESH_QUALITY is requested in VOLUME_OUTPUT and set the config boolean accordingly. ---*/
  Wrt_MeshQuality = false;
  for (unsigned short iField = 0; iField < nVolumeOutput; iField++) {
    if(VolumeOutput[iField].find("MESH_QUALITY") != string::npos) {
      Wrt_MeshQuality = true;
    }
  }

  /*--- Check if MULTIGRID is requested in VOLUME_OUTPUT and set the config boolean accordingly. ---*/
  Wrt_MultiGrid = false;
  for (unsigned short iField = 0; iField < nVolumeOutput; iField++) {
    if(VolumeOutput[iField].find("MULTIGRID") != string::npos) {
      Wrt_MultiGrid = true;
    }
  }

  if (Kind_Solver == MAIN_SOLVER::NAVIER_STOKES && Kind_Turb_Model != TURB_MODEL::NONE){
    SU2_MPI::Error("KIND_TURB_MODEL must be NONE if SOLVER= NAVIER_STOKES", CURRENT_FUNCTION);
  }
  if (Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES && Kind_Turb_Model != TURB_MODEL::NONE){
    SU2_MPI::Error("KIND_TURB_MODEL must be NONE if SOLVER= INC_NAVIER_STOKES", CURRENT_FUNCTION);
  }
  if (Kind_Solver == MAIN_SOLVER::RANS && Kind_Turb_Model == TURB_MODEL::NONE){
    SU2_MPI::Error("A turbulence model must be specified with KIND_TURB_MODEL if SOLVER= RANS", CURRENT_FUNCTION);
  }
  if (Kind_Solver == MAIN_SOLVER::INC_RANS && Kind_Turb_Model == TURB_MODEL::NONE){
    SU2_MPI::Error("A turbulence model must be specified with KIND_TURB_MODEL if SOLVER= INC_RANS", CURRENT_FUNCTION);
  }
  if (Kind_Turb_Model == TURB_MODEL::NONE && Kind_Trans_Model != TURB_TRANS_MODEL::NONE) {
    SU2_MPI::Error("KIND_TURB_MODEL cannot be NONE to use a transition model", CURRENT_FUNCTION);
  }
  switch (Kind_Solver) {
    case MAIN_SOLVER::EULER:
    case MAIN_SOLVER::INC_EULER:
    case MAIN_SOLVER::FEM_EULER:
    case MAIN_SOLVER::NEMO_EULER:
      if (nMarker_HeatFlux + nMarker_Isothermal + nMarker_HeatTransfer +
          nMarker_Smoluchowski_Maxwell + nMarker_CHTInterface > 0) {
        SU2_MPI::Error("Euler solvers are only compatible with slip walls (MARKER_EULER)", CURRENT_FUNCTION);
      }
      break;
    default:
      break;
  }

  /*--- Postprocess SST_OPTIONS into structure. ---*/
  if (Kind_Turb_Model == TURB_MODEL::SST) {
    sstParsedOptions = ParseSSTOptions(SST_Options, nSST_Options, rank);
  } else if (Kind_Turb_Model == TURB_MODEL::SA) {
    saParsedOptions = ParseSAOptions(SA_Options, nSA_Options, rank);
  }

  if (Kind_Solver == MAIN_SOLVER::INC_RANS && sstParsedOptions.compSarkar){
    SU2_MPI::Error("COMPRESSIBILITY-SARKAR only supported for SOLVER= RANS", CURRENT_FUNCTION);
  }

  if (Kind_Solver == MAIN_SOLVER::INC_RANS && sstParsedOptions.compWilcox){
    SU2_MPI::Error("COMPRESSIBILITY-WILCOX only supported for SOLVER= RANS", CURRENT_FUNCTION);
  }

  /*--- Postprocess LM_OPTIONS into structure. ---*/
  if (Kind_Trans_Model == TURB_TRANS_MODEL::LM) {
    lmParsedOptions = ParseLMOptions(LM_Options, nLM_Options, rank, Kind_Turb_Model);

    /*--- Check if problem is 2D and LM2015 has been selected ---*/
    if (lmParsedOptions.LM2015 && val_nDim == 2) {
      SU2_MPI::Error("LM2015 is available only for 3D problems", CURRENT_FUNCTION);
    }
  } else if (Kind_Trans_Model == TURB_TRANS_MODEL::AFT) {
    aftParsedOptions = ParseAFTOptions(AFT_Options, nAFT_Options, rank);
  }

  /*--- Set the boolean Wall_Functions equal to true if there is a
   definition for the wall founctions ---*/

  Wall_Functions = false;
  if (nMarker_WallFunctions > 0) {
    for (iMarker = 0; iMarker < nMarker_WallFunctions; iMarker++) {
      if (Kind_WallFunctions[iMarker] != WALL_FUNCTIONS::NONE)
        Wall_Functions = true;

      if ((Kind_WallFunctions[iMarker] == WALL_FUNCTIONS::ADAPTIVE_FUNCTION) ||
          (Kind_WallFunctions[iMarker] == WALL_FUNCTIONS::SCALABLE_FUNCTION) ||
          (Kind_WallFunctions[iMarker] == WALL_FUNCTIONS::NONEQUILIBRIUM_MODEL))
        SU2_MPI::Error(string("For RANS problems, use NONE, STANDARD_WALL_FUNCTION or EQUILIBRIUM_WALL_MODEL.\n"), CURRENT_FUNCTION);

      if (Kind_WallFunctions[iMarker] == WALL_FUNCTIONS::STANDARD_FUNCTION) {
        if ((Kind_Solver != MAIN_SOLVER::RANS) && (Kind_Solver != MAIN_SOLVER::INC_RANS))
          SU2_MPI::Error(string("Wall model STANDARD_FUNCTION only available for RANS or INC_RANS.\n"), CURRENT_FUNCTION);
        if (nRough_Wall != 0)
          SU2_MPI::Error(string("Wall model STANDARD_FUNCTION and WALL_ROUGHNESS migh not be compatible. Checking required!\n"), CURRENT_FUNCTION);
      }

    }
  }

  /*--- Initialize the AoA and Sideslip variables for the incompressible
   solver. This is typically unused (often internal flows). Also fixed CL
   mode for incompressible flows is not implemented ---*/

  if (Kind_Solver == MAIN_SOLVER::INC_EULER ||
      Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES ||
      Kind_Solver == MAIN_SOLVER::INC_RANS) {

    /*--- Compute x-velocity with a safegaurd for 0.0. ---*/

    su2double Vx = 1e-10;
    if (vel_init[0] != 0.0) {
      Vx = vel_init[0];
    }

    /*--- Compute the angle-of-attack and sideslip. ---*/

    su2double alpha = 0.0, beta = 0.0;
    if (val_nDim == 2) {
      alpha = atan(vel_init[1]/Vx)*180.0/PI_NUMBER;
    } else {
      alpha = atan(vel_init[2]/Vx)*180.0/PI_NUMBER;
      beta  = atan(vel_init[1]/Vx)*180.0/PI_NUMBER;
    }

    /*--- Set alpha and beta in the config class. ---*/

    SetAoA(alpha);
    SetAoS(beta);

    if (Fixed_CL_Mode) {
      SU2_MPI::Error(string("Fixed CL mode not implemented for the incompressible solver. \n"), CURRENT_FUNCTION);
    }

    /*--- Inc CHT simulation, but energy equation of fluid is inactive. ---*/
    if (Multizone_Problem && (nMarker_CHTInterface > 0) && !Energy_Equation)
      SU2_MPI::Error(string("You probably want to set INC_ENERGY_EQUATION= YES for the fluid solver. \n"), CURRENT_FUNCTION);
  }

  /*--- Check correctness and consistency of contact resistance options. ---*/
  if (nMarker_ContactResistance > 0) {

    /*--- Set constant contact resistance across CHT interfaces if a single value is provided. ---*/
    if (nMarker_ContactResistance == 1) {
      auto val_CHTInterface = CHT_ContactResistance[0];
      delete [] CHT_ContactResistance;
      CHT_ContactResistance = new su2double[nMarker_CHTInterface];
      for (auto iCHTMarker=0u; iCHTMarker < nMarker_CHTInterface; iCHTMarker++)
        CHT_ContactResistance[iCHTMarker] = val_CHTInterface;
    }else if((nMarker_CHTInterface/2) != nMarker_ContactResistance){
      SU2_MPI::Error("Number of CHT interfaces does not match number of contact resistances.", CURRENT_FUNCTION);
    }
    for (auto iCHTMarker=0u; iCHTMarker < nMarker_ContactResistance; iCHTMarker++){
      if (CHT_ContactResistance[iCHTMarker] < 0)
        SU2_MPI::Error("Contact resistance value should be positive.", CURRENT_FUNCTION);
    }
  }

  /*--- By default, in 2D we should use TWOD_AIRFOIL (independenly from the input file) ---*/

  if (val_nDim == 2) Geo_Description = TWOD_AIRFOIL;

  /*--- Store the SU2 module that we are executing. ---*/

  Kind_SU2 = val_software;

  /*--- Set limiter for no MUSCL reconstructions ---*/

  auto SetScalarDefaults = [](bool muscl, unsigned short& KindConvScheme, UPWIND& KindUpwind, LIMITER& KindLimiter) {
    if (KindConvScheme == NO_CONVECTIVE) {
      KindConvScheme = SPACE_UPWIND;
      KindUpwind = UPWIND::SCALAR_UPWIND;
    } else if (KindConvScheme == SPACE_CENTERED) {
      SU2_MPI::Error("Centered schemes are not available for scalar transport", CURRENT_FUNCTION);
    }
    if (!muscl) KindLimiter = LIMITER::NONE;
  };
  SetScalarDefaults(MUSCL_Turb, Kind_ConvNumScheme_Turb, Kind_Upwind_Turb, Kind_SlopeLimit_Turb);
  SetScalarDefaults(MUSCL_Heat, Kind_ConvNumScheme_Heat, Kind_Upwind_Heat, Kind_SlopeLimit_Heat);
  SetScalarDefaults(MUSCL_Species, Kind_ConvNumScheme_Species, Kind_Upwind_Species, Kind_SlopeLimit_Species);

  if (MUSCL_Flow && (Kind_ConvNumScheme_Flow == SPACE_CENTERED)) {
    if (OptionIsSet("MUSCL_FLOW")) {
      SU2_MPI::Error("Centered schemes do not use MUSCL reconstruction (use MUSCL_FLOW= NO).", CURRENT_FUNCTION);
    } else {
      MUSCL_Flow = false;
    }
  }
  if (MUSCL_AdjFlow && (Kind_ConvNumScheme_AdjFlow == SPACE_CENTERED)) {
    if (OptionIsSet("MUSCL_ADJFLOW")) {
      SU2_MPI::Error("Centered schemes do not use MUSCL reconstruction (use MUSCL_ADJFLOW= NO).", CURRENT_FUNCTION);
    } else {
      MUSCL_AdjFlow = false;
    }
  }

  if (!MUSCL_Flow || (Kind_ConvNumScheme_Flow == SPACE_CENTERED)) Kind_SlopeLimit_Flow = LIMITER::NONE;
  if (!MUSCL_AdjFlow || (Kind_ConvNumScheme_AdjFlow == SPACE_CENTERED)) Kind_SlopeLimit_AdjFlow = LIMITER::NONE;
  if (!MUSCL_AdjTurb || (Kind_ConvNumScheme_AdjTurb == SPACE_CENTERED)) Kind_SlopeLimit_AdjTurb = LIMITER::NONE;

  /*--- Set the default for thrust in ActDisk ---*/

  if ((Kind_ActDisk == NET_THRUST) || (Kind_ActDisk == BC_THRUST)
      || (Kind_ActDisk == DRAG_MINUS_THRUST) || (Kind_ActDisk == MASSFLOW)
      || (Kind_ActDisk == POWER))
    ActDisk_Jump = RATIO;

  if(Marker_ActDiskBemInlet_CG && Marker_ActDiskBemInlet_Axis){
    if(nMarker_ActDiskBemInlet_CG != nMarker_ActDiskBemInlet_Axis){
      SU2_MPI::Error("Marker lists supplied to MARKER_ACTDISK_BEM_CG and MARKER_ACTDISK_BEM_AXIS must be identical.", CURRENT_FUNCTION);
    }
    for(iMarker=0; iMarker<nMarker_ActDiskBemInlet_CG; iMarker++){
      if(Marker_ActDiskBemInlet_CG[iMarker]!=Marker_ActDiskBemInlet_Axis[iMarker]){
          SU2_MPI::Error("Marker lists supplied to MARKER_ACTDISK_BEM_CG and MARKER_ACTDISK_BEM_AXIS must be identical.", CURRENT_FUNCTION);
      }
    }
  }

  if(Marker_ActDiskBemOutlet_CG && Marker_ActDiskBemOutlet_Axis){
    if(nMarker_ActDiskBemOutlet_CG != nMarker_ActDiskBemOutlet_Axis){
      SU2_MPI::Error("Marker lists supplied to MARKER_ACTDISK_BEM_CG and MARKER_ACTDISK_BEM_AXIS must be identical.", CURRENT_FUNCTION);
    }
    for(iMarker=0; iMarker<nMarker_ActDiskBemOutlet_CG; iMarker++){
      if(Marker_ActDiskBemOutlet_CG[iMarker]!=Marker_ActDiskBemOutlet_Axis[iMarker]){
          SU2_MPI::Error("Marker lists supplied to MARKER_ACTDISK_BEM_CG and MARKER_ACTDISK_BEM_AXIS must be identical.", CURRENT_FUNCTION);
      }
    }
  }

  /*--- Error-catching and automatic array adjustments for objective, marker, and weights arrays --- */

  /*--- If Kind_Obj has not been specified, these arrays need to take a default --*/

  if (Weight_ObjFunc == nullptr && Kind_ObjFunc == nullptr) {
    Kind_ObjFunc = new unsigned short[1];
    Kind_ObjFunc[0] = DRAG_COEFFICIENT;
    Weight_ObjFunc = new su2double[1];
    Weight_ObjFunc[0] = 1.0;
    nObj=1;
    nObjW=1;
  }

  /*--- Maker sure that arrays are the same length ---*/

  if (nObj>0) {
    if (nMarker_Monitoring!=nObj && Marker_Monitoring!= nullptr) {
      if (nMarker_Monitoring==1) {
        /*-- If only one marker was listed with multiple objectives, set that marker as the marker for each objective ---*/
        nMarker_Monitoring = nObj;
        string marker = Marker_Monitoring[0];
        delete[] Marker_Monitoring;
        Marker_Monitoring = new string[nMarker_Monitoring];
        for (iMarker=0; iMarker<nMarker_Monitoring; iMarker++)
          Marker_Monitoring[iMarker] = marker;
      }
      else if(nObj==1){
        /*--- If one objective and more than one marker: repeat objective over each marker, evenly weighted ---*/
        unsigned int obj = Kind_ObjFunc[0];
        su2double wt=1.0;
        delete[] Kind_ObjFunc;
        if (Weight_ObjFunc!=nullptr){
         wt = Weight_ObjFunc[0];
         delete[] Weight_ObjFunc;
        }
        Kind_ObjFunc = new short unsigned int[nMarker_Monitoring];
        Weight_ObjFunc = new su2double[nMarker_Monitoring];
        for (unsigned short iObj=0; iObj<nMarker_Monitoring; iObj++){
          Kind_ObjFunc[iObj] = obj;
          Weight_ObjFunc[iObj] = wt;
        }
        nObjW = nObj;
      }
      else if(nObj>1) {
        SU2_MPI::Error("When using more than one OBJECTIVE_FUNCTION, MARKER_MONITORING must be the same length or length 1.\n"
                       "For multiple surfaces per objective, either use one objective or list the objective multiple times.\n"
                       "For multiple objectives per marker either use one marker or list the marker multiple times.\n"
                       "Similar rules apply for multi-objective optimization using OPT_OBJECTIVE rather than OBJECTIVE_FUNCTION.",
                       CURRENT_FUNCTION);
      }
    }
  }

  /*-- Correct for case where Weight_ObjFunc has not been provided or has length < kind_objfunc---*/

  if (nObjW<nObj) {
    if (Weight_ObjFunc!= nullptr && nObjW>1) {
      SU2_MPI::Error("The option OBJECTIVE_WEIGHT must either have the same length as OBJECTIVE_FUNCTION,\n"
                     "be lenght 1, or be deleted from the config file (equal weights will be applied).", CURRENT_FUNCTION);
    }
    Weight_ObjFunc = new su2double[nObj];
    for (unsigned short iObj=0; iObj<nObj; iObj++)
      Weight_ObjFunc[iObj] = 1.0;
  }

  /*--- One final check for multi-objective with the set of objectives
   that are not counted per-surface. We will disable multi-objective here. ---*/

  if (nObj > 1) {
    unsigned short Obj_0 = Kind_ObjFunc[0];
    for (unsigned short iObj=1; iObj<nObj; iObj++){
      switch(Kind_ObjFunc[iObj]) {
        case INVERSE_DESIGN_PRESSURE:
        case INVERSE_DESIGN_HEATFLUX:
        case THRUST_COEFFICIENT:
        case TORQUE_COEFFICIENT:
        case FIGURE_OF_MERIT:
        case SURFACE_TOTAL_PRESSURE:
        case SURFACE_STATIC_PRESSURE:
        case SURFACE_STATIC_TEMPERATURE:
        case SURFACE_MASSFLOW:
        case SURFACE_UNIFORMITY:
        case SURFACE_SECONDARY:
        case SURFACE_MOM_DISTORTION:
        case SURFACE_SECOND_OVER_UNIFORM:
        case SURFACE_PRESSURE_DROP:
        case SURFACE_SPECIES_0:
        case SURFACE_SPECIES_VARIANCE:
        case CUSTOM_OBJFUNC:
          if (Kind_ObjFunc[iObj] != Obj_0) {
            SU2_MPI::Error("The following objectives can only be used for the first surface in a multi-objective \n"
                           "problem or as a single objective applied to multiple monitoring markers:\n"
                           "INVERSE_DESIGN_PRESSURE, INVERSE_DESIGN_HEATFLUX, THRUST_COEFFICIENT, TORQUE_COEFFICIENT\n"
                           "FIGURE_OF_MERIT, SURFACE_TOTAL_PRESSURE, SURFACE_STATIC_PRESSURE, SURFACE_MASSFLOW\n"
                           "SURFACE_UNIFORMITY, SURFACE_SECONDARY, SURFACE_MOM_DISTORTION, SURFACE_SECOND_OVER_UNIFORM\n"
                           "SURFACE_PRESSURE_DROP, SURFACE_STATIC_TEMPERATURE, SURFACE_SPECIES_0\n"
                           "SURFACE_SPECIES_VARIANCE, CUSTOM_OBJFUNC.\n", CURRENT_FUNCTION);
          }
          break;
        default:
          break;
      }
    }
  }

  if (nObj > 0){
    if (Kind_ObjFunc[0] == CUSTOM_OBJFUNC && CustomObjFunc.empty() && !Multizone_Problem) {
      SU2_MPI::Error("The expression for the custom objective function was not set.\n"
                    "For example, CUSTOM_OBJFUNC= LIFT/DRAG", CURRENT_FUNCTION);
    }
  }

  /*--- Check for unsteady problem ---*/

  if ((TimeMarching == TIME_MARCHING::TIME_STEPPING ||
       TimeMarching == TIME_MARCHING::DT_STEPPING_1ST ||
       TimeMarching == TIME_MARCHING::DT_STEPPING_2ND) && !Time_Domain){
    SU2_MPI::Error("TIME_DOMAIN must be set to YES if TIME_MARCHING is "
                   "TIME_STEPPING, DUAL_TIME_STEPPING-1ST_ORDER or DUAL_TIME_STEPPING-2ND_ORDER", CURRENT_FUNCTION);
  }

  if (Time_Domain){
    Delta_UnstTime = Time_Step;

    if (TimeMarching == TIME_MARCHING::TIME_STEPPING){ InnerIter = 1; }

    /*--- Set History write freq for inner and outer iteration to zero by default, so only time iterations write. ---*/
    if (!OptionIsSet("HISTORY_WRT_FREQ_INNER")) { HistoryWrtFreq[2] = 0; }
    if (!OptionIsSet("HISTORY_WRT_FREQ_OUTER")) { HistoryWrtFreq[1] = 0; }

    if (Restart == NO) {
      Restart_Iter = 0;
    } else {
      if(nTimeIter <= Restart_Iter) SU2_MPI::Error("TIME_ITER must be larger than RESTART_ITER.", CURRENT_FUNCTION);
    }

    /*--- WINDOW_START_ITER must be larger than or equal to: RESTART_ITER. Otherwise, the running average is wrong. ---*/
    if (OptionIsSet("WINDOW_START_ITER")) {
      if (StartWindowIteration < Restart_Iter) {
        SU2_MPI::Error("WINDOW_START_ITER must be larger than or equal to: RESTART_ITER!", CURRENT_FUNCTION);
      }
    } else {
      /*--- Enforced default behavior: start of the window is the first new iteration. ---*/
      if (rank == MASTER_NODE) cout << "WARNING: Setting WINDOW_START_ITER = RESTART_ITER for meaningful running average.\n";
      StartWindowIteration = Restart_Iter;
    }

    if (Time_Step <= 0.0 && Unst_CFL == 0.0){ SU2_MPI::Error("Invalid value for TIME_STEP.", CURRENT_FUNCTION); }
  } else {
    nTimeIter = 1;
    Time_Step = 0;

    /*--- Entry 0 corresponds to unsteady simulation so for steady simulation are just set to 1. ---*/
    ScreenWrtFreq[0]  = 1;
    HistoryWrtFreq[0] = 1;

    if (TimeMarching != TIME_MARCHING::HARMONIC_BALANCE) { TimeMarching = TIME_MARCHING::STEADY; }
  }

  if (Time_Domain && !GetWrt_Restart_Overwrite()){
    SU2_MPI::Error("Appending iterations to the filename (WRT_RESTART_OVERWRITE=NO) is incompatible with transient problems.", CURRENT_FUNCTION);
  }
  if (Time_Domain && !GetWrt_Surface_Overwrite()){
    SU2_MPI::Error("Appending iterations to the filename (WRT_SURFACE_OVERWRITE=NO) is incompatible with transient problems.", CURRENT_FUNCTION);
  }
  if (Time_Domain && !GetWrt_Volume_Overwrite()){
    SU2_MPI::Error("Appending iterations to the filename (WRT_VOLUME_OVERWRITE=NO) is incompatible with transient problems.", CURRENT_FUNCTION);
  }


  /*--- Ensure that Discard_InFiles is false, owerwise the gradient could be wrong ---*/

  if ((ContinuousAdjoint || DiscreteAdjoint) && Fixed_CL_Mode && !Eval_dOF_dCX)
    Discard_InFiles = false;

  /*--- Deactivate the multigrid in the adjoint problem ---*/

  if ((ContinuousAdjoint && !MG_AdjointFlow) ||
      (TimeMarching == TIME_MARCHING::TIME_STEPPING)) { nMGLevels = 0; }

  if (Kind_Solver == MAIN_SOLVER::EULER ||
      Kind_Solver == MAIN_SOLVER::NAVIER_STOKES ||
      Kind_Solver == MAIN_SOLVER::RANS ||
      Kind_Solver == MAIN_SOLVER::NEMO_EULER ||
      Kind_Solver == MAIN_SOLVER::NEMO_NAVIER_STOKES ||
      Kind_Solver == MAIN_SOLVER::FEM_EULER ||
      Kind_Solver == MAIN_SOLVER::FEM_NAVIER_STOKES ||
      Kind_Solver == MAIN_SOLVER::FEM_RANS ||
      Kind_Solver == MAIN_SOLVER::FEM_LES){
    Kind_Regime = ENUM_REGIME::COMPRESSIBLE;
  } else if (Kind_Solver == MAIN_SOLVER::INC_EULER ||
             Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES ||
             Kind_Solver == MAIN_SOLVER::INC_RANS){
    Kind_Regime = ENUM_REGIME::INCOMPRESSIBLE;
  }  else {
    Kind_Regime = ENUM_REGIME::NO_FLOW;
  }

  if ((rank == MASTER_NODE) && ContinuousAdjoint && (Ref_NonDim == DIMENSIONAL) && (Kind_SU2 == SU2_COMPONENT::SU2_CFD)) {
    cout << "WARNING: The adjoint solver should use a non-dimensional flow solution." << endl;
  }

  /*--- Initialize non-physical points/reconstructions to zero ---*/

  Nonphys_Points   = 0;
  Nonphys_Reconstr = 0;

  /*--- Set the number of external iterations to 1 for the steady state problem ---*/

  if (Kind_Solver == MAIN_SOLVER::FEM_ELASTICITY) {
    nMGLevels = 0;
    if (Kind_Struct_Solver == STRUCT_DEFORMATION::SMALL){
      MinLogResidual = log10(Linear_Solver_Error);
    }
  }

  Radiation = (Kind_Radiation != RADIATION_MODEL::NONE);

  /*--- Check for unsupported features. ---*/

  if ((Kind_Solver != MAIN_SOLVER::EULER && Kind_Solver != MAIN_SOLVER::NAVIER_STOKES && Kind_Solver != MAIN_SOLVER::RANS) && (TimeMarching == TIME_MARCHING::HARMONIC_BALANCE)){
    SU2_MPI::Error("Harmonic Balance not yet implemented for the incompressible solver.", CURRENT_FUNCTION);
  }

  /*--- Check for Fluid model consistency ---*/

  if (standard_air) {
    if (Gamma != 1.4 || Gas_Constant != 287.058) {
      Gamma = 1.4;
      Gas_Constant = 287.058;
    }
  }

/*--- Set default values for various fluid properties. ---*/

  const su2double Molecular_Weight_Default = 28.96;
  const su2double Mu_Constant_Default = (SystemMeasurements == SI) ? 1.716E-5 : (1.716E-5 / 47.88025898);
  const su2double Mu_Ref_Default = Mu_Constant_Default;
  const su2double Mu_Temperature_Ref_Default = (SystemMeasurements == SI) ? 273.15 : (273.15 * 1.8);
  const su2double Mu_S_Default = (SystemMeasurements == SI) ? 110.4 : (110.4 * 1.8);
  const su2double Specific_Heat_Cp_Default = 1004.703;
  const su2double Thermal_Conductivity_Constant_Default = (SystemMeasurements == SI) ? 2.57E-2 : (2.57E-2 * 0.577789317);
  const su2double Prandtl_Lam_Default = 0.72;
  const su2double Prandtl_Turb_Default = 0.9;
  const su2double Lewis_Number_Default = 1.0;

  auto SetDefaultIfEmpty = [](su2double*& array, unsigned short& size, const su2double& default_val) {
    if (array == nullptr) {
      array = new su2double[1];
      array[0] = default_val;
      size = 1;
    }
  };

  SetDefaultIfEmpty(Molecular_Weight, nMolecular_Weight, Molecular_Weight_Default);
  SetDefaultIfEmpty(Specific_Heat_Cp, nSpecific_Heat_Cp, Specific_Heat_Cp_Default);
  SetDefaultIfEmpty(Mu_Constant, nMu_Constant, Mu_Constant_Default);
  if (Mu_Ref == nullptr && Mu_Temperature_Ref == nullptr && Mu_S == nullptr) {
    SetDefaultIfEmpty(Mu_Ref, nMu_Ref, Mu_Ref_Default);
    SetDefaultIfEmpty(Mu_Temperature_Ref, nMu_Temperature_Ref, Mu_Temperature_Ref_Default);
    SetDefaultIfEmpty(Mu_S, nMu_S, Mu_S_Default);
  }
  SetDefaultIfEmpty(Thermal_Conductivity_Constant, nThermal_Conductivity_Constant,
                    Thermal_Conductivity_Constant_Default);
  SetDefaultIfEmpty(Prandtl_Lam, nPrandtl_Lam, Prandtl_Lam_Default);
  SetDefaultIfEmpty(Prandtl_Turb, nPrandtl_Turb, Prandtl_Turb_Default);
  SetDefaultIfEmpty(Constant_Lewis_Number, nConstant_Lewis_Number, Lewis_Number_Default);

  Variable_Density = ((Kind_DensityModel == INC_DENSITYMODEL::VARIABLE) || (Kind_DensityModel == INC_DENSITYMODEL::FLAMELET));

  /*--- Check whether inputs for FLUID_MIXTURE are correctly specified. ---*/

    if (Kind_FluidModel == FLUID_MIXTURE) {
      /*--- Check whether the number of entries of each specified fluid property equals the number of transported scalar
       equations solved + 1. nMolecular_Weight and nSpecific_Heat_Cp are used because it is required for the fluid mixing models.
       * Cp is required in case of MIXTURE_FLUID_MODEL because the energy equation needs to be active.--- */
      if (nMolecular_Weight != nSpecies_Init + 1 || nSpecific_Heat_Cp != nSpecies_Init + 1) {
        SU2_MPI::Error(
            "The use of FLUID_MIXTURE requires the number of entries for MOLECULAR_WEIGHT and SPECIFIC_HEAT_CP,\n"
            "to be equal to the number of entries of SPECIES_INIT + 1",
            CURRENT_FUNCTION);
      }
      /*--- Check whether the density model used is correct, in the case of FLUID_MIXTURE the density model must be
       VARIABLE. Otherwise, if the density model is CONSTANT, the scalars will not have influence the mixture density
       and it will remain constant through the complete domain. --- */
      if (Kind_DensityModel != INC_DENSITYMODEL::VARIABLE) {
        SU2_MPI::Error("The use of FLUID_MIXTURE requires the INC_DENSITY_MODEL option to be VARIABLE",
                       CURRENT_FUNCTION);
      }
      /*--- Check whether the Kind scalar model used is correct, in the case of FLUID_MIXTURE the kind scalar model must
       be SPECIES_TRANSPORT. Otherwise, if the scalar model is NONE, the species transport equations will not be solved.
       --- */
      if (Kind_Species_Model != SPECIES_MODEL::SPECIES_TRANSPORT) {
        SU2_MPI::Error("The use of FLUID_MIXTURE requires the KIND_SCALAR_MODEL option to be SPECIES_TRANSPORT",
                       CURRENT_FUNCTION);
      }

      switch (Kind_ViscosityModel) {
        case VISCOSITYMODEL::CONSTANT:
          if (nMu_Constant != nSpecies_Init + 1) {
            SU2_MPI::Error(
                "The use of FLUID_MIXTURE requires the number of entries for MU_CONSTANT,\n"
                "to be equal to the number of entries of SPECIES_INIT + 1",
                CURRENT_FUNCTION);
          }
          break;
        case VISCOSITYMODEL::SUTHERLAND:
          if ((nMu_Ref != nSpecies_Init + 1) || (nMu_Temperature_Ref != nSpecies_Init + 1) ||
              (nMu_S != nSpecies_Init + 1)) {
            SU2_MPI::Error(
                "The use of FLUID_MIXTURE requires the number of entries for MU_REF, MU_T_REF and "
                "SUTHERLAND_CONSTANT,\n"
                "to be equal to the number of entries of SPECIES_INIT + 1",
                CURRENT_FUNCTION);
          }
          break;
        default:
          if (nSpecies_Init + 1 != 1) SU2_MPI::Error("Fluid mixture: viscosity model not available.", CURRENT_FUNCTION);
          break;
      }

      switch (Kind_ConductivityModel) {
        case CONDUCTIVITYMODEL::CONSTANT:
          if ((Kind_ConductivityModel_Turb == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) &&
              (Kind_Turb_Model != TURB_MODEL::NONE)) {
            if ((nThermal_Conductivity_Constant != nSpecies_Init + 1) || (nPrandtl_Turb != nSpecies_Init + 1)) {
              SU2_MPI::Error(
                  "The use of FLUID_MIXTURE requires the number of entries for THERMAL_CONDUCTIVITY_CONSTANT and "
                  "PRANDTL_TURB,\n"
                  "to be equal to the number of entries of SPECIES_INIT + 1",
                  CURRENT_FUNCTION);
            }
          } else {
            if (nThermal_Conductivity_Constant != nSpecies_Init + 1) {
              SU2_MPI::Error(
                  "The use of FLUID_MIXTURE requires the number of entries for THERMAL_CONDUCTIVITY_CONSTANT,\n"
                  "to be equal to the number of entries of SPECIES_INIT + 1",
                  CURRENT_FUNCTION);
            }
          }
          break;
        case CONDUCTIVITYMODEL::CONSTANT_PRANDTL:
          if (Kind_ConductivityModel_Turb == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
            if ((nPrandtl_Lam != nSpecies_Init + 1) || (nPrandtl_Turb != nSpecies_Init + 1)) {
              SU2_MPI::Error(
                  "The use of FLUID_MIXTURE requires the number of entries for PRANDTL_LAM and PRANDTL_TURB,\n"
                  "to be equal to the number of entries of SPECIES_INIT + 1",
                  CURRENT_FUNCTION);
            }
          } else {
            if (nPrandtl_Lam != nSpecies_Init + 1) {
              SU2_MPI::Error(
                  "The use of FLUID_MIXTURE requires the number of entries for PRANDTL_LAM,\n"
                  "to be equal to the number of entries of SPECIES_INIT + 1",
                  CURRENT_FUNCTION);
            }
          }
          break;
        default:
          if (nSpecies_Init + 1 != 1) SU2_MPI::Error("Conductivity model not available.", CURRENT_FUNCTION);
          break;
      }
    }


    if (Kind_Species_Model == SPECIES_MODEL::FLAMELET) {

      if (Kind_FluidModel != FLUID_FLAMELET) {
        SU2_MPI::Error("The use of SCALAR_MODEL= FLAMELET requires the FLUID_MODEL option to be FLUID_FLAMELET",
                       CURRENT_FUNCTION);
      }

      if (!Variable_Density) {
        SU2_MPI::Error("The use of FLUID_FLAMELET requires the INC_DENSITY_MODEL option to be VARIABLE or FLAMELET",
                       CURRENT_FUNCTION);
      }

      if (Kind_ConductivityModel != CONDUCTIVITYMODEL::FLAMELET) {
        SU2_MPI::Error("The use of FLUID_FLAMELET requires the CONDUCTIVITY_MODEL option to be FLAMELET",
                       CURRENT_FUNCTION);
      }

      if (Kind_Diffusivity_Model != DIFFUSIVITYMODEL::FLAMELET) {
        SU2_MPI::Error("The use of FLUID_FLAMELET requires the DIFFUSIVITY_MODEL option to be FLAMELET",
                       CURRENT_FUNCTION);
      }

      if (Kind_ViscosityModel != VISCOSITYMODEL::FLAMELET) {
        SU2_MPI::Error("The use of FLUID_FLAMELET requires the VISCOSITY_MODEL option to be FLAMELET",
                       CURRENT_FUNCTION);
      }

      if (Weakly_Coupled_Heat) {
        SU2_MPI::Error("The use of FLUID_FLAMELET is incompatible with WEAKLY_COUPLED_HEAT in the same zone.",
                       CURRENT_FUNCTION);
      }
    }

    /*--- Check for Measurement System ---*/

    if (SystemMeasurements == US && !standard_air) {
      SU2_MPI::Error("Only STANDARD_AIR fluid model can be used with US Measurement System", CURRENT_FUNCTION);
    }

    /* --- Check for NEMO compatibility issues ---*/
    if (Kind_FluidModel == SU2_NONEQ && (Kind_TransCoeffModel != TRANSCOEFFMODEL::WILKE && Kind_TransCoeffModel != TRANSCOEFFMODEL::SUTHERLAND && Kind_TransCoeffModel != TRANSCOEFFMODEL::GUPTAYOS) ) {
      SU2_MPI::Error("Transport model not available for NEMO solver using SU2TCLIB. Please use the WILKE, SUTHERLAND or GUPTAYOS transport model instead.", CURRENT_FUNCTION);
    }

    if (Kind_Solver == MAIN_SOLVER::NEMO_NAVIER_STOKES) {
      if (Kind_FluidModel == SU2_NONEQ && GasModel == "AIR-7" && Kind_TransCoeffModel != TRANSCOEFFMODEL::GUPTAYOS) {
        SU2_MPI::Error("Only Gupta-Yos transport model available for ionized flows using SU2TCLIB.", CURRENT_FUNCTION);
      }
    }

    if (Kind_FluidModel == MUTATIONPP && (Kind_TransCoeffModel == TRANSCOEFFMODEL::SUTHERLAND)) {
      SU2_MPI::Error("Transport model not available for NEMO solver using MUTATIONPP. Please use the WILKE, GUPTAYOS, or CHAPMANN_ENSKOG transport model instead.",
                     CURRENT_FUNCTION);
    }

    if (Kind_FluidModel == SU2_NONEQ && GasModel == "AIR-7" && nWall_Catalytic != 0) {
      SU2_MPI::Error("Catalytic wall recombination is not yet available for ionized flows in SU2_NEMO.", CURRENT_FUNCTION);
    }

    if (!ideal_gas && !nemo) {
      if (Kind_Upwind_Flow != UPWIND::ROE && Kind_Upwind_Flow != UPWIND::HLLC && Kind_Centered_Flow != CENTERED::JST) {
        SU2_MPI::Error("Only ROE Upwind, HLLC Upwind scheme, and JST scheme can be used for Non-Ideal Compressible Fluids", CURRENT_FUNCTION);
      }
    }

    if (GetBoolTurbomachinery()) {
      nBlades = new su2double[nZone];
      FreeStreamTurboNormal = new su2double[3];
    }

    /*--- Check if Giles are used with turbo markers ---*/

    if (nMarker_Giles > 0 && !GetBoolTurbomachinery()) {
      SU2_MPI::Error("Giles Boundary conditions can only be used with turbomachinery markers", CURRENT_FUNCTION);
    }

    /*--- Check if turbomachinery performance kind is specified with turbo markers ---*/
    if (GetBoolTurbomachinery() && !(nTurboMachineryKind/nZone == 1)){
      SU2_MPI::Error("Insufficient TURBO_PERF_KIND options specified with turbomachinery markers", CURRENT_FUNCTION);
    }

    /*--- Check for Boundary condition available for NICFD ---*/

    if ((!ideal_gas) && (!noneq_gas)) {
      if (nMarker_Inlet != 0) {
        SU2_MPI::Error(
            "Riemann Boundary conditions or Giles must be used for inlet and outlet with Not Ideal Compressible "
            "Fluids ",
            CURRENT_FUNCTION);
      }
      if (nMarker_Outlet != 0) {
        SU2_MPI::Error("Riemann Boundary conditions or Giles must be used outlet with Not Ideal Compressible Fluids ",
                       CURRENT_FUNCTION);
      }

      if (nMarker_FarField != 0) {
        SU2_MPI::Error("Riemann Boundary conditions or Giles must be used outlet with Not Ideal Compressible Fluids ",
                       CURRENT_FUNCTION);
      }
    }

    /*--- Check for Boundary condition available for NICF ---*/

    if (ideal_gas && (Kind_Solver != MAIN_SOLVER::INC_EULER && Kind_Solver != MAIN_SOLVER::INC_NAVIER_STOKES &&
                      Kind_Solver != MAIN_SOLVER::INC_RANS)) {
      if (SystemMeasurements == US && standard_air) {
        if (Kind_ViscosityModel != VISCOSITYMODEL::SUTHERLAND) {
          SU2_MPI::Error("Only SUTHERLAND viscosity model can be used with US Measurement", CURRENT_FUNCTION);
        }
      }
      if (Kind_ConductivityModel != CONDUCTIVITYMODEL::CONSTANT_PRANDTL) {
        SU2_MPI::Error("Only CONSTANT_PRANDTL thermal conductivity model can be used with STANDARD_AIR and IDEAL_GAS",
                       CURRENT_FUNCTION);
      }
    }
    /*--- Check for Boundary condition option agreement ---*/
  if (Kind_InitOption == REYNOLDS){
    if ((Kind_Solver == MAIN_SOLVER::NAVIER_STOKES || Kind_Solver == MAIN_SOLVER::RANS) && Reynolds <=0){
      SU2_MPI::Error("Reynolds number required for NAVIER_STOKES and RANS !!", CURRENT_FUNCTION);
    }
  }

  if (nKind_SurfaceMovement != nMarker_Moving) {
    SU2_MPI::Error("Number of SURFACE_MOVEMENT must match number of MARKER_MOVING", CURRENT_FUNCTION);
  }

  if (TimeMarching == TIME_MARCHING::TIME_STEPPING){
    nIter      = 1;
    nInnerIter  = 1;
  }

  if (!Multizone_Problem){
    ScreenWrtFreq[1]  = 0;
    HistoryWrtFreq[1] = 0;
    if (!Time_Domain){
      /*--- If not running multizone or unsteady, INNER_ITER and ITER are interchangeable,
       * but precedence will be given to INNER_ITER if both options are present. ---*/
      if (!OptionIsSet("INNER_ITER")){
        nInnerIter = nIter;
      }
    }
  }


  if ((Multizone_Problem || Time_Domain) && OptionIsSet("ITER")){
    SU2_MPI::Error("ITER must not be used when running multizone and/or unsteady problems.\n"
                   "Use TIME_ITER, OUTER_ITER or INNER_ITER to specify number of time iterations,\n"
                   "outer iterations or inner iterations, respectively.", CURRENT_FUNCTION);
  }

  /*--- If we're solving a purely steady problem with no prescribed grid
   movement (both rotating frame and moving walls can be steady), make sure that
   there is no grid motion ---*/

  if (GetGrid_Movement()){
    if ((Kind_SU2 == SU2_COMPONENT::SU2_CFD || Kind_SU2 == SU2_COMPONENT::SU2_SOL) &&
        (TimeMarching == TIME_MARCHING::STEADY && !Time_Domain)){

      if((Kind_GridMovement != ROTATING_FRAME) &&
         (Kind_GridMovement != STEADY_TRANSLATION) &&
         (Kind_GridMovement != NONE)){
        SU2_MPI::Error("Unsupported kind of grid movement for steady state problems.", CURRENT_FUNCTION);
      }
      for (iMarker = 0; iMarker < nMarker_Moving; iMarker++){
        if (Kind_SurfaceMovement[iMarker] != MOVING_WALL){
          SU2_MPI::Error("Unsupported kind of surface movement for steady state problems.", CURRENT_FUNCTION);
        }
      }
    }
  }

  /*--- The Line Search should be applied only in the deformation stage. ---*/

  if (Kind_SU2 != SU2_COMPONENT::SU2_DEF) {
    Opt_RelaxFactor = 1.0;
  }

  /*--- If it is not specified, set the mesh motion mach number
   equal to the freestream value. ---*/

  if (GetDynamic_Grid() && Mach_Motion == 0.0)
    Mach_Motion = Mach;

  /*--- Set the boolean flag if we are in a rotating frame (source term). ---*/

  Rotating_Frame = (Kind_GridMovement == ROTATING_FRAME);

  /*--- In case the grid movement parameters have not been declared in the
   config file, set them equal to zero for safety. Also check to make sure
   that for each option, a value has been declared for each moving marker. ---*/

  if (nMarker_Moving > 0){
    if (nMarkerMotion_Origin == 0){
      nMarkerMotion_Origin = 3*nMarker_Moving;
      MarkerMotion_Origin = new su2double[nMarkerMotion_Origin] ();
    }
    if (nMarkerMotion_Origin/3 != nMarker_Moving){
      SU2_MPI::Error("Number of SURFACE_MOTION_ORIGIN must be three times the number of MARKER_MOVING, (x,y,z) per marker.", CURRENT_FUNCTION);
    }
    if (nMarkerTranslation == 0){
      nMarkerTranslation = 3*nMarker_Moving;
      MarkerTranslation_Rate = new su2double[nMarkerTranslation] ();
    }
    if (nMarkerTranslation/3 != nMarker_Moving){
      SU2_MPI::Error("Number of SURFACE_TRANSLATION_RATE must be three times the number of MARKER_MOVING, (x,y,z) per marker.", CURRENT_FUNCTION);
    }
    if (nMarkerRotation_Rate == 0){
      nMarkerRotation_Rate = 3*nMarker_Moving;
      MarkerRotation_Rate = new su2double[nMarkerRotation_Rate] ();
    }
    if (nMarkerRotation_Rate/3 != nMarker_Moving){
      SU2_MPI::Error("Number of SURFACE_ROTATION_RATE must be three times the number of MARKER_MOVING, (x,y,z) per marker.", CURRENT_FUNCTION);
    }
    if (nMarkerPlunging_Ampl == 0){
      nMarkerPlunging_Ampl = 3*nMarker_Moving;
      MarkerPlunging_Ampl = new su2double[nMarkerPlunging_Ampl] ();
    }
    if (nMarkerPlunging_Ampl/3 != nMarker_Moving){
      SU2_MPI::Error("Number of SURFACE_PLUNGING_AMPL must be three times the number of MARKER_MOVING, (x,y,z) per marker.", CURRENT_FUNCTION);
    }
    if (nMarkerPlunging_Omega == 0){
      nMarkerPlunging_Omega = 3*nMarker_Moving;
      MarkerPlunging_Omega = new su2double[nMarkerPlunging_Omega] ();
    }
    if (nMarkerPlunging_Omega/3 != nMarker_Moving){
      SU2_MPI::Error("Number of SURFACE_PLUNGING_OMEGA must be three times the number of MARKER_MOVING, (x,y,z) per marker.", CURRENT_FUNCTION);
    }
    if (nMarkerPitching_Ampl == 0){
      nMarkerPitching_Ampl = 3*nMarker_Moving;
      MarkerPitching_Ampl = new su2double[nMarkerPitching_Ampl] ();
    }
    if (nMarkerPitching_Ampl/3 != nMarker_Moving){
      SU2_MPI::Error("Number of SURFACE_PITCHING_AMPL must be three times the number of MARKER_MOVING, (x,y,z) per marker.", CURRENT_FUNCTION);
    }
    if (nMarkerPitching_Omega == 0){
      nMarkerPitching_Omega = 3*nMarker_Moving;
      MarkerPitching_Omega = new su2double[nMarkerPitching_Omega] ();
    }
    if (nMarkerPitching_Omega/3 != nMarker_Moving){
      SU2_MPI::Error("Number of SURFACE_PITCHING_OMEGA must be three times the number of MARKER_MOVING, (x,y,z) per marker.", CURRENT_FUNCTION);
    }
    if (nMarkerPitching_Phase == 0){
      nMarkerPitching_Phase = 3*nMarker_Moving;
      MarkerPitching_Phase = new su2double[nMarkerPitching_Phase] ();
    }
    if (nMarkerPitching_Phase/3 != nMarker_Moving){
      SU2_MPI::Error("Number of SURFACE_PITCHING_PHASE must be three times the number of MARKER_MOVING, (x,y,z) per marker.", CURRENT_FUNCTION);
    }

    if (nMoveMotion_Origin == 0){
      nMoveMotion_Origin = nMarker_Moving;
      MoveMotion_Origin = new unsigned short[nMoveMotion_Origin];
      for (iMarker = 0; iMarker < nMarker_Moving; iMarker++){
        MoveMotion_Origin[iMarker] = NO;
      }
    }
    if (nMoveMotion_Origin != nMarker_Moving){
      SU2_MPI::Error("Number of MOVE_MOTION_ORIGIN must match number of MARKER_MOVING.", CURRENT_FUNCTION);
    }
  }

  /*-- Setting Harmonic Balance period from the config file */

  if (TimeMarching == TIME_MARCHING::HARMONIC_BALANCE) {
    HarmonicBalance_Period = GetHarmonicBalance_Period();
    if (HarmonicBalance_Period < 0)  {
      SU2_MPI::Error("Not a valid value for time period!!", CURRENT_FUNCTION);
    }
    /* Initialize the Harmonic balance Frequency pointer */
    if (Omega_HB == nullptr) {
      Omega_HB = new su2double[nOmega_HB];
      for (unsigned short iZone = 0; iZone < nOmega_HB; iZone++ )
        Omega_HB[iZone] = 0.0;
  } else {
      if (nOmega_HB != nTimeInstances) {
        SU2_MPI::Error("Length of omega_HB  must match the number TIME_INSTANCES!!" , CURRENT_FUNCTION);
      }
    }
  }

  /*--- Force number of span-wise section to 1 if 2D case ---*/
  if(val_nDim ==2){
    nSpanWiseSections_User=1;
    Kind_SpanWise= EQUISPACED;
  }

  /*--- Set number of TurboPerformance markers ---*/
  if(nMarker_Turbomachinery > 0){
    if(nMarker_Turbomachinery > 1){
      nMarker_TurboPerformance = nMarker_Turbomachinery + SU2_TYPE::Int(nMarker_Turbomachinery/2) + 1;
    }else{
      nMarker_TurboPerformance = nMarker_Turbomachinery;
    }
  } else {
    nMarker_TurboPerformance = 0;
    nSpanWiseSections =1;
  }

  /*--- Set number of TurboPerformance markers ---*/
  if(nMarker_Turbomachinery != 0){
    nSpan_iZones = new unsigned short[nZone];
  }

  /*--- Set number of TurboPerformance markers ---*/
  if(GetGrid_Movement() && RampRotatingFrame && !DiscreteAdjoint){
    FinalRotation_Rate_Z = Rotation_Rate[2];
    if(abs(FinalRotation_Rate_Z) > 0.0){
      Rotation_Rate[2] = rampRotFrame_coeff[0];
    }
  }

  if(RampOutletPressure && !DiscreteAdjoint){
    for (iMarker = 0; iMarker < nMarker_Giles; iMarker++){
      if (Kind_Data_Giles[iMarker] == STATIC_PRESSURE || Kind_Data_Giles[iMarker] == STATIC_PRESSURE_1D || Kind_Data_Giles[iMarker] == RADIAL_EQUILIBRIUM ){
        FinalOutletPressure = Giles_Var1[iMarker];
        Giles_Var1[iMarker] = rampOutPres_coeff[0];
      }
    }
    for (iMarker = 0; iMarker < nMarker_Riemann; iMarker++){
      if (Kind_Data_Riemann[iMarker] == STATIC_PRESSURE || Kind_Data_Riemann[iMarker] == RADIAL_EQUILIBRIUM){
        FinalOutletPressure = Riemann_Var1[iMarker];
        Riemann_Var1[iMarker] = rampOutPres_coeff[0];
      }
    }
  }

  /*--- Check on extra Relaxation factor for Giles---*/
  if(extrarelfac[1] > 0.5){
    extrarelfac[1] = 0.5;
  }
    /*--- Use the various rigid-motion input frequencies to determine the period to be used with harmonic balance cases.
     There are THREE types of motion to consider, namely: rotation, pitching, and plunging.
     The largest period of motion is the one to be used for harmonic balance  calculations. ---*/

  /*if (Unsteady_Simulation == HARMONIC_BALANCE) {
    if (!(GetGrid_Movement())) {
      // No grid movement - Time period from config file //
      HarmonicBalance_Period = GetHarmonicBalance_Period();
    }

    else {
      unsigned short N_MOTION_TYPES = 3;
      su2double *periods;
      periods = new su2double[N_MOTION_TYPES];

      //--- rotation: ---//

      su2double Omega_mag_rot = sqrt(pow(Rotation_Rate_X[ZONE_0],2)+pow(Rotation_Rate_Y[ZONE_0],2)+pow(Rotation_Rate_Z[ZONE_0],2));
      if (Omega_mag_rot > 0)
          periods[0] = 2*PI_NUMBER/Omega_mag_rot;
      else
          periods[0] = 0.0;

      //--- pitching: ---//

      su2double Omega_mag_pitch = sqrt(pow(Pitching_Omega_X[ZONE_0],2)+pow(Pitching_Omega_Y[ZONE_0],2)+pow(Pitching_Omega_Z[ZONE_0],2));
      if (Omega_mag_pitch > 0)
          periods[1] = 2*PI_NUMBER/Omega_mag_pitch;
      else
          periods[1] = 0.0;

      //--- plunging: ---//

      su2double Omega_mag_plunge = sqrt(pow(Plunging_Omega_X[ZONE_0],2)+pow(Plunging_Omega_Y[ZONE_0],2)+pow(Plunging_Omega_Z[ZONE_0],2));
      if (Omega_mag_plunge > 0)
          periods[2] = 2*PI_NUMBER/Omega_mag_plunge;
      else
          periods[2] = 0.0;

      //--- determine which period is largest ---//

      unsigned short iVar;
      HarmonicBalance_Period = 0.0;
      for (iVar = 0; iVar < N_MOTION_TYPES; iVar++) {
          if (periods[iVar] > HarmonicBalance_Period)
              HarmonicBalance_Period = periods[iVar];
      }

      delete periods;
    }

  }*/


  /*--- In case the moment origin coordinates have not been declared in the
   config file, set them equal to zero for safety. Also check to make sure
   that for each marker, a value has been declared for the moment origin.
   Unless only one value was specified, then set this value for all the markers
   being monitored. ---*/


  if ((nRefOriginMoment_X != nRefOriginMoment_Y) || (nRefOriginMoment_X != nRefOriginMoment_Z) ) {
    SU2_MPI::Error("ERROR: Length of REF_ORIGIN_MOMENT_X, REF_ORIGIN_MOMENT_Y and REF_ORIGIN_MOMENT_Z must be the same!!", CURRENT_FUNCTION);
  }

  if (RefOriginMoment_X == nullptr) {
    RefOriginMoment_X = new su2double[nMarker_Monitoring];
    for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
      RefOriginMoment_X[iMarker] = 0.0;
  } else {
    if (nRefOriginMoment_X == 1) {

      su2double aux_RefOriginMoment_X = RefOriginMoment_X[0];
      delete [] RefOriginMoment_X;
      RefOriginMoment_X = new su2double[nMarker_Monitoring];
      nRefOriginMoment_X = nMarker_Monitoring;

      for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
        RefOriginMoment_X[iMarker] = aux_RefOriginMoment_X;
    }
    else if (nRefOriginMoment_X != nMarker_Monitoring) {
      SU2_MPI::Error("ERROR: Length of REF_ORIGIN_MOMENT_X must match number of Monitoring Markers!!", CURRENT_FUNCTION);
    }
  }

  if (RefOriginMoment_Y == nullptr) {
    RefOriginMoment_Y = new su2double[nMarker_Monitoring];
    for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
      RefOriginMoment_Y[iMarker] = 0.0;
  } else {
    if (nRefOriginMoment_Y == 1) {

      su2double aux_RefOriginMoment_Y = RefOriginMoment_Y[0];
      delete [] RefOriginMoment_Y;
      RefOriginMoment_Y = new su2double[nMarker_Monitoring];
      nRefOriginMoment_Y = nMarker_Monitoring;

      for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
        RefOriginMoment_Y[iMarker] = aux_RefOriginMoment_Y;
    }
    else if (nRefOriginMoment_Y != nMarker_Monitoring) {
      SU2_MPI::Error("ERROR: Length of REF_ORIGIN_MOMENT_Y must match number of Monitoring Markers!!", CURRENT_FUNCTION);
    }
  }

  if (RefOriginMoment_Z == nullptr) {
    RefOriginMoment_Z = new su2double[nMarker_Monitoring];
    for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
      RefOriginMoment_Z[iMarker] = 0.0;
  } else {
    if (nRefOriginMoment_Z == 1) {

      su2double aux_RefOriginMoment_Z = RefOriginMoment_Z[0];
      delete [] RefOriginMoment_Z;
      RefOriginMoment_Z = new su2double[nMarker_Monitoring];
      nRefOriginMoment_Z = nMarker_Monitoring;

      for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ )
        RefOriginMoment_Z[iMarker] = aux_RefOriginMoment_Z;
    }
    else if (nRefOriginMoment_Z != nMarker_Monitoring) {
      SU2_MPI::Error("ERROR: Length of REF_ORIGIN_MOMENT_Z must match number of Monitoring Markers!!", CURRENT_FUNCTION);
    }
  }

  /*--- Set the boolean flag if we are carrying out an aeroelastic simulation. ---*/

  Aeroelastic_Simulation = GetGrid_Movement() && (GetSurface_Movement(AEROELASTIC) || GetSurface_Movement(AEROELASTIC_RIGID_MOTION));

  /*--- Initializing the size for the solutions of the Aeroelastic problem. ---*/


  if (GetGrid_Movement() && Aeroelastic_Simulation) {
    Aeroelastic_np1.resize(nMarker_Monitoring);
    Aeroelastic_n.resize(nMarker_Monitoring);
    Aeroelastic_n1.resize(nMarker_Monitoring);
    for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++) {
      Aeroelastic_np1[iMarker].resize(2);
      Aeroelastic_n[iMarker].resize(2);
      Aeroelastic_n1[iMarker].resize(2);
      for (int i =0; i<2; i++) {
        Aeroelastic_np1[iMarker][i].resize(2);
        Aeroelastic_n[iMarker][i].resize(2);
        Aeroelastic_n1[iMarker][i].resize(2);
        for (int j=0; j<2; j++) {
          Aeroelastic_np1[iMarker][i][j] = 0.0;
          Aeroelastic_n[iMarker][i][j] = 0.0;
          Aeroelastic_n1[iMarker][i][j] = 0.0;
        }
      }
    }
  }

  /*--- Allocate memory for the plunge and pitch and initialized them to zero ---*/

  if (GetGrid_Movement() && Aeroelastic_Simulation) {
    Aeroelastic_pitch = new su2double[nMarker_Monitoring];
    Aeroelastic_plunge = new su2double[nMarker_Monitoring];
    for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++ ) {
      Aeroelastic_pitch[iMarker] = 0.0;
      Aeroelastic_plunge[iMarker] = 0.0;
    }
  }

  FinestMesh = MESH_0;
  if (MGCycle == FULLMG_CYCLE) FinestMesh = nMGLevels;

  if ((Kind_Solver == MAIN_SOLVER::NAVIER_STOKES) &&
      (Kind_Turb_Model != TURB_MODEL::NONE))
    Kind_Solver = MAIN_SOLVER::RANS;

  if ((Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES) &&
      (Kind_Turb_Model != TURB_MODEL::NONE))
    Kind_Solver = MAIN_SOLVER::INC_RANS;

  if (Kind_Solver == MAIN_SOLVER::EULER ||
      Kind_Solver == MAIN_SOLVER::INC_EULER ||
      Kind_Solver == MAIN_SOLVER::NEMO_EULER ||
      Kind_Solver == MAIN_SOLVER::FEM_EULER)
    Kind_Turb_Model = TURB_MODEL::NONE;

  Kappa_2nd_Flow = jst_coeff[0];
  Kappa_4th_Flow = jst_coeff[1];
  Kappa_2nd_AdjFlow = jst_adj_coeff[0];
  Kappa_4th_AdjFlow = jst_adj_coeff[1];

  /*--- Make the MG_PreSmooth, MG_PostSmooth, and MG_CorrecSmooth
   arrays consistent with nMGLevels ---*/

  auto * tmp_smooth = new unsigned short[nMGLevels+1];

  if ((nMG_PreSmooth != nMGLevels+1) && (nMG_PreSmooth != 0)) {
    if (nMG_PreSmooth > nMGLevels+1) {

      /*--- Truncate by removing unnecessary elements at the end ---*/

      for (unsigned int i = 0; i <= nMGLevels; i++)
        tmp_smooth[i] = MG_PreSmooth[i];
      delete [] MG_PreSmooth;
      MG_PreSmooth=nullptr;
    }
    else {

      /*--- Add additional elements equal to last element ---*/

      for (unsigned int i = 0; i < nMG_PreSmooth; i++)
        tmp_smooth[i] = MG_PreSmooth[i];
      for (unsigned int i = nMG_PreSmooth; i <= nMGLevels; i++)
        tmp_smooth[i] = MG_PreSmooth[nMG_PreSmooth-1];
      delete [] MG_PreSmooth;
      MG_PreSmooth=nullptr;
    }

    nMG_PreSmooth = nMGLevels+1;
    MG_PreSmooth = new unsigned short[nMG_PreSmooth];
    for (unsigned int i = 0; i < nMG_PreSmooth; i++)
      MG_PreSmooth[i] = tmp_smooth[i];
  }
  if ((nMGLevels != 0) && (nMG_PreSmooth == 0)) {
    delete [] MG_PreSmooth;
    nMG_PreSmooth = nMGLevels+1;
    MG_PreSmooth = new unsigned short[nMG_PreSmooth];
    for (unsigned int i = 0; i < nMG_PreSmooth; i++)
      MG_PreSmooth[i] = i+1;
  }

  if ((nMG_PostSmooth != nMGLevels+1) && (nMG_PostSmooth != 0)) {
    if (nMG_PostSmooth > nMGLevels+1) {

      /*--- Truncate by removing unnecessary elements at the end ---*/

      for (unsigned int i = 0; i <= nMGLevels; i++)
        tmp_smooth[i] = MG_PostSmooth[i];
      delete [] MG_PostSmooth;
      MG_PostSmooth=nullptr;
    }
    else {

      /*--- Add additional elements equal to last element ---*/

      for (unsigned int i = 0; i < nMG_PostSmooth; i++)
        tmp_smooth[i] = MG_PostSmooth[i];
      for (unsigned int i = nMG_PostSmooth; i <= nMGLevels; i++)
        tmp_smooth[i] = MG_PostSmooth[nMG_PostSmooth-1];
      delete [] MG_PostSmooth;
      MG_PostSmooth=nullptr;
    }

    nMG_PostSmooth = nMGLevels+1;
    MG_PostSmooth = new unsigned short[nMG_PostSmooth];
    for (unsigned int i = 0; i < nMG_PostSmooth; i++)
      MG_PostSmooth[i] = tmp_smooth[i];

  }

  if ((nMGLevels != 0) && (nMG_PostSmooth == 0)) {
    delete [] MG_PostSmooth;
    nMG_PostSmooth = nMGLevels+1;
    MG_PostSmooth = new unsigned short[nMG_PostSmooth];
    for (unsigned int i = 0; i < nMG_PostSmooth; i++)
      MG_PostSmooth[i] = 0;
  }

  if ((nMG_CorrecSmooth != nMGLevels+1) && (nMG_CorrecSmooth != 0)) {
    if (nMG_CorrecSmooth > nMGLevels+1) {

      /*--- Truncate by removing unnecessary elements at the end ---*/

      for (unsigned int i = 0; i <= nMGLevels; i++)
        tmp_smooth[i] = MG_CorrecSmooth[i];
      delete [] MG_CorrecSmooth;
      MG_CorrecSmooth = nullptr;
    }
    else {

      /*--- Add additional elements equal to last element ---*/

      for (unsigned int i = 0; i < nMG_CorrecSmooth; i++)
        tmp_smooth[i] = MG_CorrecSmooth[i];
      for (unsigned int i = nMG_CorrecSmooth; i <= nMGLevels; i++)
        tmp_smooth[i] = MG_CorrecSmooth[nMG_CorrecSmooth-1];
      delete [] MG_CorrecSmooth;
      MG_CorrecSmooth = nullptr;
    }
    nMG_CorrecSmooth = nMGLevels+1;
    MG_CorrecSmooth = new unsigned short[nMG_CorrecSmooth];
    for (unsigned int i = 0; i < nMG_CorrecSmooth; i++)
      MG_CorrecSmooth[i] = tmp_smooth[i];
  }

  if ((nMGLevels != 0) && (nMG_CorrecSmooth == 0)) {
    delete [] MG_CorrecSmooth;
    nMG_CorrecSmooth = nMGLevels+1;
    MG_CorrecSmooth = new unsigned short[nMG_CorrecSmooth];
    for (unsigned int i = 0; i < nMG_CorrecSmooth; i++)
      MG_CorrecSmooth[i] = 0;
  }

  /*--- Override MG Smooth parameters ---*/

  if (nMG_PreSmooth != 0) MG_PreSmooth[MESH_0] = 1;
  if (nMG_PostSmooth != 0) {
    MG_PostSmooth[MESH_0] = 0;
    MG_PostSmooth[nMGLevels] = 0;
  }
  if (nMG_CorrecSmooth != 0) MG_CorrecSmooth[nMGLevels] = 0;

  if (Restart) MGCycle = V_CYCLE;

  if (ContinuousAdjoint) {
    if (Kind_Solver == MAIN_SOLVER::EULER) Kind_Solver = MAIN_SOLVER::ADJ_EULER;
    if (Kind_Solver == MAIN_SOLVER::NAVIER_STOKES) Kind_Solver = MAIN_SOLVER::ADJ_NAVIER_STOKES;
    if (Kind_Solver == MAIN_SOLVER::RANS) Kind_Solver = MAIN_SOLVER::ADJ_RANS;
  }

  nCFL = nMGLevels+1;
  CFL = new su2double[nCFL];
  CFL[0] = CFLFineGrid;

  /*--- Handle optional CFL adapt parameter values ---*/

  if (nCFL_AdaptParam < default_cfl_adapt.size()) {
    auto newParam = new su2double [default_cfl_adapt.size()];
    for (iCFL = 0; iCFL < default_cfl_adapt.size(); ++iCFL) {
      if (iCFL < nCFL_AdaptParam) newParam[iCFL] = CFL_AdaptParam[iCFL];
      else newParam[iCFL] = default_cfl_adapt[iCFL];
    }
    swap(newParam, CFL_AdaptParam);
    delete [] newParam;
    nCFL_AdaptParam = default_cfl_adapt.size();
  }

  /*--- Evaluate when the Cl should be evaluated ---*/

  Iter_Fixed_CM        = SU2_TYPE::Int(nInnerIter / (su2double(Update_iH)+1));
  Iter_Fixed_NetThrust = SU2_TYPE::Int(nInnerIter / (su2double(Update_BCThrust)+1));

  /*--- Setting relaxation factor and CFL for the adjoint runs ---*/

  if (ContinuousAdjoint) {
    CFL[0] = CFL[0] * CFLRedCoeff_AdjFlow;
    CFL_AdaptParam[2] *= CFLRedCoeff_AdjFlow;
    CFL_AdaptParam[3] *= CFLRedCoeff_AdjFlow;
    Iter_Fixed_CM = SU2_TYPE::Int(su2double (Iter_Fixed_CM) / CFLRedCoeff_AdjFlow);
    Iter_Fixed_NetThrust = SU2_TYPE::Int(su2double (Iter_Fixed_NetThrust) / CFLRedCoeff_AdjFlow);
  }

  if ((DiscreteAdjoint) && (Inconsistent_Disc)) {
    Kind_ConvNumScheme_Flow = Kind_ConvNumScheme_AdjFlow;
    Kind_Centered_Flow = Kind_Centered_AdjFlow;
    Kind_Upwind_Flow = Kind_Upwind_AdjFlow;
    Kappa_2nd_Flow = jst_adj_coeff[0];
    Kappa_4th_Flow = jst_adj_coeff[1];
  }

  if (Update_AoA_Iter_Limit == 0 && Fixed_CL_Mode) {
    SU2_MPI::Error("ERROR: Please specify non-zero UPDATE_AOA_ITER_LIMIT.", CURRENT_FUNCTION);
  }
  if (Iter_Fixed_CM == 0) { Iter_Fixed_CM = nInnerIter+1; Update_iH = 0; }
  if (Iter_Fixed_NetThrust == 0) { Iter_Fixed_NetThrust = nInnerIter+1; Update_BCThrust = 0; }

  for (iCFL = 1; iCFL < nCFL; iCFL++)
    CFL[iCFL] = CFL[iCFL-1];

  if (nRKStep == 0) {
    nRKStep = 1;
    RK_Alpha_Step = new su2double[1]; RK_Alpha_Step[0] = 1.0;
  }

  /* Check if the byte alignment of the matrix multiplications is a
     multiple of 64. */
  if( byteAlignmentMatMul%64 ) {
    SU2_MPI::Error("ALIGNED_BYTES_MATMUL must be a multiple of 64.", CURRENT_FUNCTION);
  }

  /* Determine the value of sizeMatMulPadding, which is the matrix size in
     the vectorization direction when padding is applied to have optimal
     performance in the matrix multiplications. */
  sizeMatMulPadding = byteAlignmentMatMul/sizeof(passivedouble);

  /* Correct the number of time levels for time accurate local time
     stepping, if needed.  */
  if (nLevels_TimeAccurateLTS == 0)  nLevels_TimeAccurateLTS =  1;
  if (nLevels_TimeAccurateLTS  > 15) nLevels_TimeAccurateLTS = 15;

  /* Check that no time accurate local time stepping is specified for time
     integration schemes other than ADER. */
  if (Kind_TimeIntScheme_FEM_Flow != ADER_DG && nLevels_TimeAccurateLTS != 1) {

    if (rank==MASTER_NODE) {
      cout << endl << "WARNING: "
           << nLevels_TimeAccurateLTS << " levels specified for time accurate local time stepping." << endl
           << "Time accurate local time stepping is only possible for ADER, hence this option is not used." << endl
           << endl;
    }

    nLevels_TimeAccurateLTS = 1;
  }

  if (Kind_TimeIntScheme_FEM_Flow == ADER_DG) {

    TimeMarching = TIME_MARCHING::TIME_STEPPING;  // Only time stepping for ADER.

    /* If time accurate local time stepping is used, make sure that an unsteady
       CFL is specified. If not, terminate. */
    if (nLevels_TimeAccurateLTS != 1) {
      if(Unst_CFL == 0.0)
        SU2_MPI::Error("ERROR: Unsteady CFL not specified for time accurate local time stepping.",
                       CURRENT_FUNCTION);
    }

    /* Determine the location of the ADER time DOFs, which are the Gauss-Legendre
       integration points corresponding to the number of time DOFs. */
    vector<passivedouble> GLPoints(nTimeDOFsADER_DG), GLWeights(nTimeDOFsADER_DG);
    CGaussJacobiQuadrature GaussJacobi;
    GaussJacobi.GetQuadraturePoints(0.0, 0.0, -1.0, 1.0, GLPoints, GLWeights);

    TimeDOFsADER_DG = new su2double[nTimeDOFsADER_DG];
    for(unsigned short i=0; i<nTimeDOFsADER_DG; ++i)
      TimeDOFsADER_DG[i] = GLPoints[i];

    /* Determine the number of integration points in time, their locations
       on the interval [-1..1] and their integration weights. */
    unsigned short orderExact = ceil(Quadrature_Factor_Time_ADER_DG*(nTimeDOFsADER_DG-1));
    nTimeIntegrationADER_DG = orderExact/2 + 1;
    nTimeIntegrationADER_DG = max(nTimeIntegrationADER_DG, nTimeDOFsADER_DG);
    GLPoints.resize(nTimeIntegrationADER_DG);
    GLWeights.resize(nTimeIntegrationADER_DG);
    GaussJacobi.GetQuadraturePoints(0.0, 0.0, -1.0, 1.0, GLPoints, GLWeights);

    TimeIntegrationADER_DG    = new su2double[nTimeIntegrationADER_DG];
    WeightsIntegrationADER_DG = new su2double[nTimeIntegrationADER_DG];
    for(unsigned short i=0; i<nTimeIntegrationADER_DG; ++i) {
      TimeIntegrationADER_DG[i]    = GLPoints[i];
      WeightsIntegrationADER_DG[i] = GLWeights[i];
    }
  }

  if(Kind_TimeIntScheme_Turb != EULER_IMPLICIT &&
     Kind_TimeIntScheme_Turb != EULER_EXPLICIT){
    SU2_MPI::Error("Only TIME_DISCRE_TURB = EULER_IMPLICIT, EULER_EXPLICIT have been implemented.", CURRENT_FUNCTION);
  }

  if (nIntCoeffs == 0) {
    nIntCoeffs = 2;
    Int_Coeffs = new su2double[2]; Int_Coeffs[0] = 0.25; Int_Coeffs[1] = 0.5;
  }

  if (nElasticityMod == 0) {
    nElasticityMod = 1;
    ElasticityMod = new su2double[1]; ElasticityMod[0] = 2E11;
  }

  if (nPoissonRatio == 0) {
    nPoissonRatio = 1;
    PoissonRatio = new su2double[1]; PoissonRatio[0] = 0.30;
  }

  if (nMaterialDensity == 0) {
    nMaterialDensity = 1;
    MaterialDensity = new su2double[1]; MaterialDensity[0] = 7854;
  }

  if (nMaterialThermalExpansion == 0) {
    nMaterialThermalExpansion = 1;
    MaterialThermalExpansion = new su2double[1]();
  }

  if (nElasticityMod != nPoissonRatio || nElasticityMod != nMaterialDensity ||
      nElasticityMod != nMaterialThermalExpansion) {
    SU2_MPI::Error("ELASTICITY_MODULUS, POISSON_RATIO, MATERIAL_DENSITY, and THERMAL_EXPANSION_COEFF need "
                   "to have the same number of entries (the number of materials).", CURRENT_FUNCTION);
  }

  if (nElectric_Constant == 0) {
    nElectric_Constant = 1;
    Electric_Constant = new su2double[1]; Electric_Constant[0] = 0.0;
  }

  if (nElectric_Field == 0) {
    nElectric_Field = 1;
    Electric_Field_Mod = new su2double[1]; Electric_Field_Mod[0] = 0.0;
  }

  if (nDim_RefNode == 0) {
    nDim_RefNode = 3;
    RefNode_Displacement = new su2double[3];
    RefNode_Displacement[0] = 0.0; RefNode_Displacement[1] = 0.0; RefNode_Displacement[2] = 0.0;
  }

  if (nDim_Electric_Field == 0) {
    nDim_Electric_Field = 2;
    Electric_Field_Dir = new su2double[2]; Electric_Field_Dir[0] = 0.0;  Electric_Field_Dir[1] = 1.0;
  }

  if ((Kind_SU2 == SU2_COMPONENT::SU2_CFD) && (Kind_Solver == MAIN_SOLVER::NONE)) {
    SU2_MPI::Error("PHYSICAL_PROBLEM must be set in the configuration file", CURRENT_FUNCTION);
  }

  /*--- Set a flag for viscous simulations ---*/

  Viscous = (( Kind_Solver == MAIN_SOLVER::NAVIER_STOKES          ) ||
             ( Kind_Solver == MAIN_SOLVER::NEMO_NAVIER_STOKES     ) ||
             ( Kind_Solver == MAIN_SOLVER::ADJ_NAVIER_STOKES      ) ||
             ( Kind_Solver == MAIN_SOLVER::RANS                   ) ||
             ( Kind_Solver == MAIN_SOLVER::ADJ_RANS               ) ||
             ( Kind_Solver == MAIN_SOLVER::FEM_NAVIER_STOKES      ) ||
             ( Kind_Solver == MAIN_SOLVER::FEM_RANS               ) ||
             ( Kind_Solver == MAIN_SOLVER::FEM_LES                ) ||
             ( Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES      ) ||
             ( Kind_Solver == MAIN_SOLVER::INC_RANS               ) );

  /*--- To avoid boundary intersections, let's add a small constant to the planes. ---*/

  if (Geo_Description == NACELLE) {
    for (unsigned short iSections = 0; iSections < nLocationStations; iSections++) {
      if (LocationStations[iSections] == 0) LocationStations[iSections] = 1E-6;
      if (LocationStations[iSections] == 360) LocationStations[iSections] = 359.999999;
    }
  }
  else {
    for (unsigned short iSections = 0; iSections < nLocationStations; iSections++) {
      LocationStations[iSections] += EPS;
    }
    geo_loc[0] += EPS;
    geo_loc[1] += EPS;
  }

  /*--- Length based parameter for slope limiters uses a default value of
   1.0m ---*/

  RefElemLength = 1.0;
  if (SystemMeasurements == US) RefElemLength /= 0.3048;

  /*--- Re-scale the length based parameters. The US system uses feet,
   but SU2 assumes that the grid is in inches ---*/

  if ((SystemMeasurements == US) && (Kind_SU2 == SU2_COMPONENT::SU2_CFD)) {

    for (iMarker = 0; iMarker < nMarker_Monitoring; iMarker++) {
      RefOriginMoment_X[iMarker] = RefOriginMoment_X[iMarker]/12.0;
      RefOriginMoment_Y[iMarker] = RefOriginMoment_Y[iMarker]/12.0;
      RefOriginMoment_Z[iMarker] = RefOriginMoment_Z[iMarker]/12.0;
    }

    for (iMarker = 0; iMarker < nMarker_Moving; iMarker++){
      for (unsigned short iDim = 0; iDim < 3; iDim++){
        MarkerMotion_Origin[3*iMarker+iDim] /= 12.0;
      }
    }

    RefLength = RefLength/12.0;

    if ((val_nDim == 2) && (!Axisymmetric)) RefArea = RefArea/12.0;
    else RefArea = RefArea/144.0;
    Length_Reynolds = Length_Reynolds/12.0;
    Highlite_Area = Highlite_Area/144.0;
    SemiSpan = SemiSpan/12.0;

    ea_lim[0] /= 12.0;
    ea_lim[1] /= 12.0;
    ea_lim[2] /= 12.0;

    if (Geo_Description != NACELLE) {
      for (unsigned short iSections = 0; iSections < nLocationStations; iSections++) {
        LocationStations[iSections] = LocationStations[iSections]/12.0;
      }
      geo_loc[0] /= 12.0;
      geo_loc[1] /= 12.0;
    }

    for (int i=0; i<7; ++i) eng_cyl[i] /= 12.0;
  }

  if(Turb_Fixed_Values && !OptionIsSet("TURB_FIXED_VALUES_DOMAIN")){
    SU2_MPI::Error("TURB_FIXED_VALUES activated, but no domain set with TURB_FIXED_VALUES_DOMAIN.", CURRENT_FUNCTION);
  }

  /*--- Check for constant lift mode. Initialize the update flag for
   the AoA with each iteration to false  ---*/

  if (Fixed_CL_Mode) Update_AoA = false;

  if (DirectDiff != NO_DERIVATIVE) {
#ifndef CODI_FORWARD_TYPE
    if (Kind_SU2 == SU2_COMPONENT::SU2_CFD) {
      SU2_MPI::Error("SU2_CFD: Config option DIRECT_DIFF= YES requires AD support.\n"
                     "Please use SU2_CFD_DIRECTDIFF (meson.py ... -Denable-directdiff=true ...).",
                     CURRENT_FUNCTION);
    }
#endif
    /*--- Initialize the derivative values ---*/
    switch (DirectDiff) {
      case D_MACH:
        SU2_TYPE::SetDerivative(Mach, 1.0);
        break;
      case D_AOA:
        SU2_TYPE::SetDerivative(AoA, 1.0);
        break;
      case D_SIDESLIP:
        SU2_TYPE::SetDerivative(AoS, 1.0);
        break;
      case D_REYNOLDS:
        SU2_TYPE::SetDerivative(Reynolds, 1.0);
        break;
      case D_TURB2LAM:
       SU2_TYPE::SetDerivative(TurbIntensityAndViscRatioFreeStream[1], 1.0);
        break;
      default:
        /*--- All other cases are handled in the specific solver ---*/
        break;
      }
  }

#if defined CODI_REVERSE_TYPE
  AD_Mode = YES;

  AD::PreaccEnabled = AD_Preaccumulation;

#else
  if (AD_Mode == YES) {
    SU2_MPI::Error("Config option AUTO_DIFF= YES requires AD support.\n"
                   "Please use SU2_???_AD (meson.py ... -Denable-autodiff=true ...).",
                   CURRENT_FUNCTION);
  }
#endif

  delete [] tmp_smooth;

  /*--- Make sure that implicit time integration is disabled
        for the FEM fluid solver (numerics). ---*/
  if ((Kind_Solver == MAIN_SOLVER::FEM_EULER)         ||
      (Kind_Solver == MAIN_SOLVER::FEM_NAVIER_STOKES) ||
      (Kind_Solver == MAIN_SOLVER::FEM_RANS)          ||
      (Kind_Solver == MAIN_SOLVER::FEM_LES)) {
     Kind_TimeIntScheme_Flow = Kind_TimeIntScheme_FEM_Flow;
  }

  /*--- Set up the time stepping / unsteady CFL options. ---*/
  if ((TimeMarching == TIME_MARCHING::TIME_STEPPING) && (Unst_CFL != 0.0)) {
    for (iCFL = 0; iCFL < nCFL; iCFL++)
      CFL[iCFL] = Unst_CFL;
  }


  /*--- If it is a fixed mode problem, then we will add Iter_dCL_dAlpha iterations to
    evaluate the derivatives with respect to a change in the AoA and CL ---*/

  if (!ContinuousAdjoint & !DiscreteAdjoint) {
    if (Fixed_CL_Mode) nInnerIter += Iter_dCL_dAlpha;
  }

  /* --- Set Finite Difference mode to false by default --- */

  Finite_Difference_Mode = false;

  /*--- If there are not design variables defined in the file ---*/

  if (nDV == 0) {
    nDV = 1;
    Design_Variable = new unsigned short [nDV];
    Design_Variable[0] = NO_DEFORMATION;
  }

  /*--- Checks for incompressible flow problems. ---*/

  if (Kind_Solver == MAIN_SOLVER::INC_EULER) {
    /*--- Force inviscid problems to use constant density and disable energy. ---*/
    if (Kind_DensityModel != INC_DENSITYMODEL::CONSTANT || Energy_Equation) {
      SU2_MPI::Error("Inviscid incompressible problems must be constant density (no energy eqn.).\n Use DENSITY_MODEL= CONSTANT and ENERGY_EQUATION= NO.", CURRENT_FUNCTION);
    }
  }

  /*--- Default values should recover original incompressible behavior (for old config files). ---*/

  if (Kind_Solver == MAIN_SOLVER::INC_EULER || Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES || Kind_Solver == MAIN_SOLVER::INC_RANS) {
    if ((Kind_DensityModel == INC_DENSITYMODEL::CONSTANT) || (Kind_DensityModel == INC_DENSITYMODEL::BOUSSINESQ))
      Kind_FluidModel = CONSTANT_DENSITY;
  }

  if ((Kind_DensityModel != INC_DENSITYMODEL::CONSTANT) && (Kind_Species_Model==SPECIES_MODEL::NONE)) Energy_Equation = true;
  /*--- For the flamelet combustion model, energy equation is a passive field, we lookup T and write it to the field ---*/
  if (Kind_Species_Model == SPECIES_MODEL::FLAMELET ) Energy_Equation = false;

  if (Kind_DensityModel == INC_DENSITYMODEL::BOUSSINESQ) {
    Energy_Equation = true;
    if (Body_Force) {
      SU2_MPI::Error("Body force and Boussinesq source terms are not currently compatible.", CURRENT_FUNCTION);
    }
  }

  if (Kind_DensityModel == INC_DENSITYMODEL::VARIABLE) {
    if (Kind_FluidModel != INC_IDEAL_GAS && Kind_FluidModel != INC_IDEAL_GAS_POLY && Kind_FluidModel != FLUID_MIXTURE && Kind_FluidModel != FLUID_FLAMELET) {
      SU2_MPI::Error("Variable density incompressible solver limited to ideal gases.\n Check the fluid model options (use INC_IDEAL_GAS, INC_IDEAL_GAS_POLY).", CURRENT_FUNCTION);
    }
  }

  if (Kind_Solver != MAIN_SOLVER::INC_EULER && Kind_Solver != MAIN_SOLVER::INC_NAVIER_STOKES && Kind_Solver != MAIN_SOLVER::INC_RANS) {
    if ((Kind_FluidModel == CONSTANT_DENSITY) || (Kind_FluidModel == INC_IDEAL_GAS) || (Kind_FluidModel == INC_IDEAL_GAS_POLY)) {
      SU2_MPI::Error("Fluid model not compatible with compressible flows.\n CONSTANT_DENSITY/INC_IDEAL_GAS/INC_IDEAL_GAS_POLY are for incompressible only.", CURRENT_FUNCTION);
    }
  }

  if (Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES || Kind_Solver == MAIN_SOLVER::INC_RANS) {
    if (Kind_ViscosityModel == VISCOSITYMODEL::SUTHERLAND) {
      if ((Kind_FluidModel != INC_IDEAL_GAS) && (Kind_FluidModel != INC_IDEAL_GAS_POLY) && (Kind_FluidModel != FLUID_MIXTURE)) {
        SU2_MPI::Error("Sutherland's law only valid for ideal gases in incompressible flows.\n Must use VISCOSITY_MODEL=CONSTANT_VISCOSITY and set viscosity with\n MU_CONSTANT, or use DENSITY_MODEL= VARIABLE with FLUID_MODEL= INC_IDEAL_GAS or INC_IDEAL_GAS_POLY for VISCOSITY_MODEL=SUTHERLAND.\n NOTE: FREESTREAM_VISCOSITY is no longer used for incompressible flows!", CURRENT_FUNCTION);
      }
    }
  }

  /*--- Vorticity confinement feature currently not supported for incompressible or non-equilibrium model or axisymmetric flows. ---*/

  if ((Kind_Solver == MAIN_SOLVER::INC_EULER
    || Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES
    || Kind_Solver == MAIN_SOLVER::INC_RANS
    || Kind_Solver == MAIN_SOLVER::NEMO_EULER
    || Kind_Solver == MAIN_SOLVER::NEMO_NAVIER_STOKES
    || Axisymmetric)
    && VorticityConfinement) {
    SU2_MPI::Error("Vorticity confinement feature currently not supported for incompressible or non-equilibrium model or axisymmetric flows.", CURRENT_FUNCTION);
  }

  /*--- Actuator disk BEM method for propellers feature currently not supported for incompressible or non-equilibrium model or axisymmetric flows. ---*/

  if ((Kind_Solver == MAIN_SOLVER::INC_EULER
    || Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES
    || Kind_Solver == MAIN_SOLVER::INC_RANS
    || Kind_Solver == MAIN_SOLVER::NEMO_EULER
    || Kind_Solver == MAIN_SOLVER::NEMO_NAVIER_STOKES
    || Axisymmetric)
    && ActDisk_DoubleSurface) {
    SU2_MPI::Error("Actuator disk BEM method for propellers feature currently not supported for incompressible or non-equilibrium model or axisymmetric flows.", CURRENT_FUNCTION);
  }

  /*--- Check the coefficients for the polynomial models. ---*/

  if (Kind_Solver != MAIN_SOLVER::INC_EULER && Kind_Solver != MAIN_SOLVER::INC_NAVIER_STOKES && Kind_Solver != MAIN_SOLVER::INC_RANS) {
    if ((Kind_ViscosityModel == VISCOSITYMODEL::POLYNOMIAL) || (Kind_ConductivityModel == CONDUCTIVITYMODEL::POLYNOMIAL) || (Kind_FluidModel == INC_IDEAL_GAS_POLY)) {
      SU2_MPI::Error("POLYNOMIAL_VISCOSITY and POLYNOMIAL_CONDUCTIVITY are for incompressible only currently.", CURRENT_FUNCTION);
    }
  }

  /*--- Data-driven fluid model is currently only supported for compressible flow problems. ---*/
  if ((Kind_Solver == MAIN_SOLVER::INC_EULER || Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES || Kind_Solver == MAIN_SOLVER::INC_RANS) && (Kind_FluidModel == DATADRIVEN_FLUID)) {
    SU2_MPI::Error("Data-driven fluid model can only be used for compressible flows.", CURRENT_FUNCTION);
  }

  if ((Kind_Solver == MAIN_SOLVER::INC_EULER || Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES || Kind_Solver == MAIN_SOLVER::INC_RANS) && (Kind_FluidModel == INC_IDEAL_GAS_POLY)) {
    su2double sum = 0.0;
    for (unsigned short iVar = 0; iVar < N_POLY_COEFFS; iVar++) {
      sum += GetCp_PolyCoeff(iVar);
    }
    if ((N_POLY_COEFFS < 1) || (sum == 0.0))
      SU2_MPI::Error(string("CP_POLYCOEFFS not set for fluid model INC_IDEAL_GAS_POLY. \n"), CURRENT_FUNCTION);
  }

  if (((Kind_Solver == MAIN_SOLVER::INC_EULER || Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES || Kind_Solver == MAIN_SOLVER::INC_RANS)) && (Kind_ViscosityModel == VISCOSITYMODEL::POLYNOMIAL)) {
    su2double sum = 0.0;
    for (unsigned short iVar = 0; iVar < N_POLY_COEFFS; iVar++) {
      sum += GetMu_PolyCoeff(iVar);
    }
    if ((N_POLY_COEFFS < 1) || (sum == 0.0))
      SU2_MPI::Error(string("MU_POLYCOEFFS not set for viscosity model POLYNOMIAL_VISCOSITY. \n"), CURRENT_FUNCTION);
  }

  if ((Kind_Solver == MAIN_SOLVER::INC_EULER || Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES || Kind_Solver == MAIN_SOLVER::INC_RANS) && (Kind_ConductivityModel == CONDUCTIVITYMODEL::POLYNOMIAL)) {
    su2double sum = 0.0;
    for (unsigned short iVar = 0; iVar < N_POLY_COEFFS; iVar++) {
      sum += GetKt_PolyCoeff(iVar);
    }
    if ((N_POLY_COEFFS < 1) || (sum == 0.0))
      SU2_MPI::Error(string("KT_POLYCOEFFS not set for conductivity model POLYNOMIAL_CONDUCTIVITY. \n"), CURRENT_FUNCTION);
  }

  /*--- Incompressible solver currently limited to SI units. ---*/

  if ((Kind_Solver == MAIN_SOLVER::INC_EULER || Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES || Kind_Solver == MAIN_SOLVER::INC_RANS) && (SystemMeasurements == US)) {
    SU2_MPI::Error("Must use SI units for incompressible solver.", CURRENT_FUNCTION);
  }

  /*--- Check that the non-dim type is valid. ---*/

  if ((Kind_Solver == MAIN_SOLVER::INC_EULER || Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES || Kind_Solver == MAIN_SOLVER::INC_RANS)) {
    if ((Ref_Inc_NonDim != INITIAL_VALUES) && (Ref_Inc_NonDim != REFERENCE_VALUES) && (Ref_Inc_NonDim != DIMENSIONAL)) {
      SU2_MPI::Error("Incompressible non-dim. scheme invalid.\n Must use INITIAL_VALUES, REFERENCE_VALUES, or DIMENSIONAL.", CURRENT_FUNCTION);
    }
  }

  /*--- Check that the incompressible inlets are correctly specified. ---*/

  if ((Kind_Solver == MAIN_SOLVER::INC_EULER || Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES || Kind_Solver == MAIN_SOLVER::INC_RANS) && (nMarker_Inlet != 0)) {
    if (nMarker_Inlet != nInc_Inlet) {
      SU2_MPI::Error("Inlet types for incompressible problem improperly specified.\n Use INC_INLET_TYPE= VELOCITY_INLET or PRESSURE_INLET.\n Must list a type for each inlet marker, including duplicates, e.g.,\n INC_INLET_TYPE= VELOCITY_INLET VELOCITY_INLET PRESSURE_INLET", CURRENT_FUNCTION);
    }
    for (unsigned short iInlet = 0; iInlet < nInc_Inlet; iInlet++){
      if ((Kind_Inc_Inlet[iInlet] != INLET_TYPE::VELOCITY_INLET) && (Kind_Inc_Inlet[iInlet] != INLET_TYPE::PRESSURE_INLET)) {
        SU2_MPI::Error("Undefined incompressible inlet type. VELOCITY_INLET or PRESSURE_INLET possible.", CURRENT_FUNCTION);
      }
    }
  }

  /*--- Check that the incompressible inlets are correctly specified. ---*/

  if ((Kind_Solver == MAIN_SOLVER::INC_EULER || Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES || Kind_Solver == MAIN_SOLVER::INC_RANS) && (nMarker_Outlet != 0)) {
    if (nMarker_Outlet != nInc_Outlet) {
      SU2_MPI::Error("Outlet types for incompressible problem improperly specified.\n Use INC_OUTLET_TYPE= PRESSURE_OUTLET or MASS_FLOW_OUTLET.\n Must list a type for each inlet marker, including duplicates, e.g.,\n INC_OUTLET_TYPE= PRESSURE_OUTLET PRESSURE_OUTLET MASS_FLOW_OUTLET", CURRENT_FUNCTION);
    }
    for (unsigned short iInlet = 0; iInlet < nInc_Outlet; iInlet++){
      if ((Kind_Inc_Outlet[iInlet] != INC_OUTLET_TYPE::PRESSURE_OUTLET) && (Kind_Inc_Outlet[iInlet] != INC_OUTLET_TYPE::MASS_FLOW_OUTLET)) {
        SU2_MPI::Error("Undefined incompressible outlet type. PRESSURE_OUTLET or MASS_FLOW_OUTLET possible.", CURRENT_FUNCTION);
      }
    }
  }

  /*--- Assert that there are two markers being analyzed if the
   pressure drop objective function is selected. ---*/

  for (unsigned short iObj = 0; iObj < nObj; iObj++) {
    if ((Kind_ObjFunc[iObj] == SURFACE_PRESSURE_DROP) && (nMarker_Analyze < 2)) {
      SU2_MPI::Error("Must list the first two markers for the pressure drop objective function.\n Expected format: MARKER_ANALYZE= (outlet_name, inlet_name, ...).", CURRENT_FUNCTION);
    }
  }

  /*--- Check feasibility for Streamwise Periodic flow ---*/
  if (Kind_Streamwise_Periodic != ENUM_STREAMWISE_PERIODIC::NONE) {
    if (Kind_Regime != ENUM_REGIME::INCOMPRESSIBLE)
      SU2_MPI::Error("Streamwise Periodic Flow currently only implemented for incompressible flow.", CURRENT_FUNCTION);
    if (Kind_Solver == MAIN_SOLVER::INC_EULER)
      SU2_MPI::Error("Streamwise Periodic Flow + Incompressible Euler: Not tested yet.", CURRENT_FUNCTION);
    if (nMarker_PerBound == 0)
      SU2_MPI::Error("A MARKER_PERIODIC pair has to be set with KIND_STREAMWISE_PERIODIC != NONE.", CURRENT_FUNCTION);
    if (Energy_Equation && Streamwise_Periodic_Temperature && nMarker_Isothermal != 0)
      SU2_MPI::Error("No MARKER_ISOTHERMAL marker allowed with STREAMWISE_PERIODIC_TEMPERATURE= YES, only MARKER_HEATFLUX & MARKER_SYM.", CURRENT_FUNCTION);
    if (Ref_Inc_NonDim != DIMENSIONAL)
      SU2_MPI::Error("Streamwise Periodicity only works with \"INC_NONDIM= DIMENSIONAL\", the nondimensionalization with source terms doesn;t work in general.", CURRENT_FUNCTION);
    if (Axisymmetric)
      SU2_MPI::Error("Streamwise Periodicity terms does not not have axisymmetric corrections.", CURRENT_FUNCTION);
    if (!Energy_Equation) Streamwise_Periodic_Temperature = false;
  } else {
    /*--- Safety measure ---*/
    Streamwise_Periodic_Temperature = false;
  }

  if (nRough_Wall > 0) {
    /*--- Validate name of the markers. ---*/
    for (iMarker = 0; iMarker < nRough_Wall; ++iMarker) {
      auto CheckMarker = [&](unsigned short nMarker, const string* markerName) {
        for (auto jMarker = 0u; jMarker < nMarker; ++jMarker) {
          if (markerName[jMarker].compare(Marker_RoughWall[iMarker]) == 0) {
            return true;
          }
        }
        return false;
      };
      if (!CheckMarker(nMarker_HeatFlux, Marker_HeatFlux) &&
          !CheckMarker(nMarker_Isothermal, Marker_Isothermal) &&
          !CheckMarker(nMarker_HeatTransfer, Marker_HeatTransfer) &&
          !CheckMarker(nMarker_CHTInterface, Marker_CHTInterface)) {
        SU2_MPI::Error("Marker " + Marker_RoughWall[iMarker] + " is not a viscous wall.", CURRENT_FUNCTION);
      }
    }
  }

  /*--- Handle default options for topology optimization ---*/

  if (topology_optimization && top_optim_nKernel==0) {
    top_optim_nKernel = 1;
    top_optim_kernels = new ENUM_FILTER_KERNEL [1];
    top_optim_kernels[0] = ENUM_FILTER_KERNEL::CONICAL_WEIGHT;
  }

  if (top_optim_nKernel != 0) {
    /*--- Set default value of kernel parameters ---*/
    if (top_optim_nKernelParams == 0) {
      top_optim_nKernelParams = top_optim_nKernel;
      top_optim_kernel_params = new su2double [top_optim_nKernel];
      for (unsigned short i=0; i<top_optim_nKernel; ++i) top_optim_kernel_params[i] = 1.0;
    }
    /*--- Broadcast the only value provided ---*/
    else if (top_optim_nKernelParams==1 && top_optim_nKernel>1) {
      su2double tmp = top_optim_kernel_params[0];
      delete [] top_optim_kernel_params;
      top_optim_nKernelParams = top_optim_nKernel;
      top_optim_kernel_params = new su2double [top_optim_nKernel];
      for (unsigned short i=0; i<top_optim_nKernel; ++i) top_optim_kernel_params[i] = tmp;
    }
    /*--- Numbers do not match ---*/
    else if (top_optim_nKernelParams != top_optim_nKernel) {
      SU2_MPI::Error("Different number of topology filter kernels and respective parameters.", CURRENT_FUNCTION);
    }

    /*--- Set default value of filter radius ---*/
    if (top_optim_nRadius == 0) {
      top_optim_nRadius = top_optim_nKernel;
      top_optim_filter_radius = new su2double [top_optim_nKernel];
      for (unsigned short i=0; i<top_optim_nKernel; ++i) top_optim_filter_radius[i] = 1.0e-6;
    }
    /*--- Broadcast the only value provided ---*/
    else if (top_optim_nRadius==1 && top_optim_nKernel>1) {
      su2double tmp = top_optim_filter_radius[0];
      delete [] top_optim_filter_radius;
      top_optim_nRadius = top_optim_nKernel;
      top_optim_filter_radius = new su2double [top_optim_nKernel];
      for (unsigned short i=0; i<top_optim_nKernel; ++i) top_optim_filter_radius[i] = tmp;
    }
    /*--- Numbers do not match ---*/
    else if (top_optim_nRadius != top_optim_nKernel) {
      SU2_MPI::Error("Different number of topology filter kernels and respective radii.", CURRENT_FUNCTION);
    }
  }

  /*--- If we are executing SU2_DOT in surface file mode, then
   force the projected surface sensitivity file to be written. ---*/

  Wrt_Projected_Sensitivity = false;
  if ((Kind_SU2 == SU2_COMPONENT::SU2_DOT) && (Design_Variable[0] == SURFACE_FILE)) {
    Wrt_Projected_Sensitivity = true;
  }

  /*--- Delay the output until exit for minimal communication mode. ---*/

  if (Comm_Level != COMM_FULL) {

    /*--- Disable the use of Comm_Level = NONE until we have properly
     implemented it. ---*/

    if (Comm_Level == COMM_NONE)
      SU2_MPI::Error("COMM_LEVEL = NONE not yet implemented.", CURRENT_FUNCTION);
  }

  /*--- Check the conductivity model. Deactivate the turbulent component
   if we are not running RANS. ---*/

  if ((Kind_Solver != MAIN_SOLVER::RANS) &&
      (Kind_Solver != MAIN_SOLVER::ADJ_RANS) &&
      (Kind_Solver != MAIN_SOLVER::DISC_ADJ_RANS) &&
      (Kind_Solver != MAIN_SOLVER::INC_RANS) &&
      (Kind_Solver != MAIN_SOLVER::DISC_ADJ_INC_RANS)){
    Kind_ConductivityModel_Turb = CONDUCTIVITYMODEL_TURB::NONE;
  }

  /* Set a default for the size of the RECTANGLE / BOX grid sizes. */

  if (nMesh_Box_Size == 0) {
    nMesh_Box_Size = 3;
    Mesh_Box_Size = new short [nMesh_Box_Size];
    Mesh_Box_Size[0] = 33;
    Mesh_Box_Size[1] = 33;
    Mesh_Box_Size[2] = 33;
  } else if (nMesh_Box_Size != 3) {
    SU2_MPI::Error("MESH_BOX_SIZE specified without 3 values.\n", CURRENT_FUNCTION);
  }

  /* Force the lowest memory preconditioner when direct solvers are used. */

  auto isPastix = [](unsigned short kindSolver) {
    return kindSolver == PASTIX_LDLT || kindSolver == PASTIX_LU;
  };

  if (isPastix(Kind_Linear_Solver)) Kind_Linear_Solver_Prec = LU_SGS;
  if (isPastix(Kind_DiscAdj_Linear_Solver)) Kind_DiscAdj_Linear_Prec = LU_SGS;
  if (isPastix(Kind_Deform_Linear_Solver)) Kind_Deform_Linear_Solver_Prec = LU_SGS;


  if (DiscreteAdjoint) {
#if !defined CODI_REVERSE_TYPE
    if (Kind_SU2 == SU2_COMPONENT::SU2_CFD) {
      SU2_MPI::Error("SU2_CFD: Config option MATH_PROBLEM= DISCRETE_ADJOINT requires AD support!\n"
                     "Please use SU2_CFD_AD (configuration/compilation is done using the preconfigure.py script).",
                     CURRENT_FUNCTION);
    }
#endif

    /*--- Use the same linear solver on the primal as the one used in the adjoint. ---*/
    Kind_Linear_Solver = Kind_DiscAdj_Linear_Solver;
    Kind_Linear_Solver_Prec = Kind_DiscAdj_Linear_Prec;

    if (Time_Domain) {

      Restart_Flow = false;

      if (Unst_AdjointIter- long(nTimeIter) < 0){
        SU2_MPI::Error("Invalid iteration number requested for unsteady adjoint.\n"
                       "Make sure EXT_ITER is larger or equal than UNST_ADJOINT_ITER.",
                       CURRENT_FUNCTION);
      }

      /*--- If the averaging interval is not set, we average over all time-steps ---*/

      if (Iter_Avg_Objective == 0.0) {
        Iter_Avg_Objective = nTimeIter;
      }

    }

    /*--- Note that this is deliberately done at the end of this routine! ---*/
    switch(Kind_Solver) {
      case MAIN_SOLVER::EULER:
        Kind_Solver = MAIN_SOLVER::DISC_ADJ_EULER;
        break;
      case MAIN_SOLVER::RANS:
        Kind_Solver = MAIN_SOLVER::DISC_ADJ_RANS;
        break;
      case MAIN_SOLVER::NAVIER_STOKES:
        Kind_Solver = MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES;
        break;
      case MAIN_SOLVER::INC_EULER:
        Kind_Solver = MAIN_SOLVER::DISC_ADJ_INC_EULER;
        break;
      case MAIN_SOLVER::INC_RANS:
        Kind_Solver = MAIN_SOLVER::DISC_ADJ_INC_RANS;
        break;
      case MAIN_SOLVER::INC_NAVIER_STOKES:
        Kind_Solver = MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES;
        break;
      case MAIN_SOLVER::FEM_EULER :
        Kind_Solver = MAIN_SOLVER::DISC_ADJ_FEM_EULER;
        break;
      case MAIN_SOLVER::FEM_RANS :
        Kind_Solver = MAIN_SOLVER::DISC_ADJ_FEM_RANS;
        break;
      case MAIN_SOLVER::FEM_NAVIER_STOKES :
        Kind_Solver = MAIN_SOLVER::DISC_ADJ_FEM_NS;
        break;
      case MAIN_SOLVER::FEM_ELASTICITY:
        Kind_Solver = MAIN_SOLVER::DISC_ADJ_FEM;
        break;
      case MAIN_SOLVER::HEAT_EQUATION:
        Kind_Solver = MAIN_SOLVER::DISC_ADJ_HEAT;
        break;
      default:
        break;
    }

    RampOutletPressure = false;
    RampRotatingFrame = false;
  }

  /* 2nd-order MUSCL is not possible for the continuous adjoint
   turbulence model. */

  if (MUSCL_AdjTurb) {
    SU2_MPI::Error("MUSCL_ADJTURB= YES not currently supported.\nPlease select MUSCL_ADJTURB= NO (first-order).",
                   CURRENT_FUNCTION);
  }

  /* Check for whether we need a second gradient method to calculate
   gradients for uwpind reconstruction. Set additional booleans to
   minimize overhead as appropriate. */

  if (MUSCL_Flow || MUSCL_Turb || MUSCL_Species || MUSCL_Heat || MUSCL_AdjFlow) {

    ReconstructionGradientRequired = true;

    if ((Kind_Gradient_Method_Recon == NO_GRADIENT) ||
        (Kind_Gradient_Method_Recon == Kind_Gradient_Method)) {

      /* The default behavior if no reconstruction gradient is specified
       is to use the same gradient as needed for the viscous/source terms
       without recomputation. If they are using the same method, then
       we also want to avoid recomputation. */

      ReconstructionGradientRequired = false;
      Kind_Gradient_Method_Recon = Kind_Gradient_Method;
    }

  } else {
    ReconstructionGradientRequired = false;
  }

  if (ReconstructionGradientRequired && GetFluidProblem() && Kind_ConvNumScheme_Flow == SPACE_CENTERED)
    SU2_MPI::Error("For centered schemes the option NUM_METHOD_GRAD_RECON should not be set.", CURRENT_FUNCTION);

  /* Simpler boolean to control allocation of least-squares memory. */

  LeastSquaresRequired = false;
  if ((Kind_Gradient_Method_Recon == LEAST_SQUARES) ||
      (Kind_Gradient_Method_Recon == WEIGHTED_LEAST_SQUARES) ||
      (Kind_Gradient_Method       == LEAST_SQUARES) ||
      (Kind_Gradient_Method       == WEIGHTED_LEAST_SQUARES)) {
    LeastSquaresRequired = true;
  }

  if (Kind_Gradient_Method == LEAST_SQUARES) {
    SU2_MPI::Error(string("LEAST_SQUARES gradient method not allowed for viscous / source terms.\n") +
                   string("Please select either WEIGHTED_LEAST_SQUARES or GREEN_GAUSS."),
                   CURRENT_FUNCTION);
  }

  /* Protect against using CFL adaption for non-flow or certain
   unsteady flow problems. */

  if (CFL_Adapt && !GetFluidProblem()) {
    SU2_MPI::Error(string("CFL adaption only available for finite-volume fluid solvers.\n") +
                   string("Please select CFL_ADAPT = NO."),
                   CURRENT_FUNCTION);
  }

  if (CFL_Adapt && (TimeMarching == TIME_MARCHING::TIME_STEPPING)) {
    SU2_MPI::Error(string("CFL adaption not available for TIME_STEPPING integration.\n") +
                   string("Please select CFL_ADAPT = NO."),
                   CURRENT_FUNCTION);
  }

  /* Protect against using incorrect CFL adaption parameters. */

  if (CFL_Adapt && (CFL_AdaptParam[0] > 1.0)) {
    SU2_MPI::Error(string("CFL adaption factor down should be less than 1.0."), CURRENT_FUNCTION);
  }

  if (CFL_Adapt && (CFL_AdaptParam[1] < 1.0)) {
    SU2_MPI::Error(string("CFL adaption factor up should be greater than 1.0."), CURRENT_FUNCTION);
  }

  if (CFL_Adapt && (CFL_AdaptParam[2] > CFL_AdaptParam[3])) {
    SU2_MPI::Error(string("CFL adaption minimum CFL is larger than the maximum CFL."), CURRENT_FUNCTION);
  }

  /*--- 0 in the config file means "disable" which can be done using a very large group. ---*/
  if (edgeColorGroupSize==0) edgeColorGroupSize = 1<<30;

  /*--- Specifying a deforming surface requires a mesh deformation solver. ---*/
  if (GetSurface_Movement(DEFORMING)) Deform_Mesh = true;

  monoatomic = GetGasModel() == "ARGON";

  /*--- Set number of Turbulence Variables. ---*/
  switch (TurbModelFamily(Kind_Turb_Model)) {
    case TURB_FAMILY::NONE:
      nTurbVar = 0; break;
    case TURB_FAMILY::SA:
      nTurbVar = 1; break;
    case TURB_FAMILY::KW:
      nTurbVar = 2; break;
  }
  /*--- Check whether the number of entries of the MARKER_INLET_TURBULENT equals the number of turbulent properties
       used for the respective turbulent model. nTurb_Properties must be equal to 1 or 2 depending on whether SA or
       SST model are used.--- */
  if (Marker_Inlet_Turb != nullptr && Kind_Turb_Model == TURB_MODEL::SST && nTurb_Properties != 2)
    SU2_MPI::Error(
        "The use of MARKER_INLET_TURBULENT requires the number of entries when SST Model is used \n"
        "to be equal to 2 : Turbulent intensity and ratio turbulent to laminar viscosity",
        CURRENT_FUNCTION);
  if (Marker_Inlet_Turb != nullptr && Kind_Turb_Model == TURB_MODEL::SA && nTurb_Properties != 1)
    SU2_MPI::Error(
        "The use of MARKER_INLET_TURBULENT requires the number of entries when SA Model is used \n"
        "to be equal to 1 : ratio turbulent to laminar viscosity",
        CURRENT_FUNCTION);

  /*--- Checks for additional species transport. ---*/
  if ((Kind_Species_Model == SPECIES_MODEL::SPECIES_TRANSPORT) || (Kind_Species_Model == SPECIES_MODEL::FLAMELET)) {
    if (Kind_Solver != MAIN_SOLVER::INC_NAVIER_STOKES &&
        Kind_Solver != MAIN_SOLVER::INC_RANS &&
        Kind_Solver != MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES &&
        Kind_Solver != MAIN_SOLVER::DISC_ADJ_INC_RANS &&
        Kind_Solver != MAIN_SOLVER::NAVIER_STOKES &&
        Kind_Solver != MAIN_SOLVER::RANS &&
        Kind_Solver != MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES &&
        Kind_Solver != MAIN_SOLVER::DISC_ADJ_RANS &&
        Kind_Solver != MAIN_SOLVER::MULTIPHYSICS)
      SU2_MPI::Error("Species transport currently only available for compressible and incompressible flow.", CURRENT_FUNCTION);

    /*--- Species specific OF currently can only handle one entry in Marker_Analyze. ---*/
    for (unsigned short iObj = 0; iObj < nObj; iObj++) {
      if ((Kind_ObjFunc[iObj] == SURFACE_SPECIES_0 ||
           Kind_ObjFunc[iObj] == SURFACE_SPECIES_VARIANCE) &&
          nMarker_Analyze > 1) {
        SU2_MPI::Error("SURFACE_SPECIES_0 and SURFACE_SPECIES_VARIANCE currently can only handle one entry to MARKER_ANALYZE.", CURRENT_FUNCTION);
      }
    }

    if(Kind_TimeIntScheme_Species != EULER_IMPLICIT &&
       Kind_TimeIntScheme_Species != EULER_EXPLICIT){
      SU2_MPI::Error("Only TIME_DISCRE_TURB = EULER_IMPLICIT, EULER_EXPLICIT have been implemented in the scalar solver.", CURRENT_FUNCTION);
    }

    /*--- If Species clipping is on, make sure bounds are given by the user. ---*/
    if (Species_Clipping)
      if (!(OptionIsSet("SPECIES_CLIPPING_MIN") && OptionIsSet("SPECIES_CLIPPING_MAX")))
        SU2_MPI::Error("SPECIES_CLIPPING= YES requires the options SPECIES_CLIPPING_MIN/MAX to set the clipping values.", CURRENT_FUNCTION);

    /*--- Make sure a Diffusivity has been set for Constant Diffusivity. ---*/
    if (Kind_Diffusivity_Model == DIFFUSIVITYMODEL::CONSTANT_DIFFUSIVITY &&
        !(OptionIsSet("DIFFUSIVITY_CONSTANT")))
      SU2_MPI::Error("A DIFFUSIVITY_CONSTANT=<value> has to be set with DIFFUSIVITY_MODEL= CONSTANT_DIFFUSIVITY.", CURRENT_FUNCTION);

    /*--- Check whether the number of entries of the constant Lewis number equals the number of transported scalar
       equations solved. nConstant_Lewis_Number is used because it is required for the diffusivity fluid mixing
       models--- */
    if (Kind_Diffusivity_Model == DIFFUSIVITYMODEL::CONSTANT_LEWIS && nConstant_Lewis_Number != nSpecies_Init + 1)
      SU2_MPI::Error(
          "The use of CONSTANT_LEWIS requires the number of entries for CONSTANT_LEWIS_NUMBER ,\n"
          "to be equal to the number of entries of SPECIES_INIT +1",
          CURRENT_FUNCTION);

    // Helper function that checks scalar variable bounds,
    auto checkScalarBounds = [&](su2double scalar, const string& name, su2double lowerBound, su2double upperBound) {
      if (scalar < lowerBound || scalar > upperBound)
        SU2_MPI::Error(string("Variable: ") + name + string(", is out of bounds."), CURRENT_FUNCTION);
    };

    /*--- Some options have to provide as many entries as there are additional species equations. ---*/
    /*--- Fill a vector with the entires and then check if each element is equal to the first one. ---*/
    std::vector<unsigned short> nSpecies_options;
    nSpecies_options.push_back(nSpecies_Init);
    if (Species_Clipping)
      nSpecies_options.insert(nSpecies_options.end(), {nSpecies_Clipping_Min, nSpecies_Clipping_Max});
    if (nMarker_Inlet_Species > 0)
      nSpecies_options.push_back(nSpecies_per_Inlet);
    // Add more options for size check here.

    /*--- nSpecies_Init is the master, but it simply checks for consistency. ---*/
    for (auto elem : nSpecies_options)
      if (nSpecies_options[0] != elem)
        SU2_MPI::Error("Make sure all species inputs have the same size.", CURRENT_FUNCTION);

    /*--- Once consistency is checked set the var that is used throughout the code. ---*/
    nSpecies = nSpecies_Init;

    /*--- Check whether some variables (or their sums) are in physical bounds. [0,1] for species related quantities. ---*/
    /*--- Note, only for species transport, not for flamelet model ---*/
    if (Kind_Species_Model == SPECIES_MODEL::SPECIES_TRANSPORT) {
      su2double Species_Init_Sum = 0.0;
      for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        checkScalarBounds(Species_Init[iSpecies], "SPECIES_INIT individual", 0.0, 1.0);
        Species_Init_Sum += Species_Init[iSpecies];
      }
      checkScalarBounds(Species_Init_Sum, "SPECIES_INIT sum", 0.0, 1.0);

      for (iMarker = 0; iMarker < nMarker_Inlet_Species; iMarker++) {
        su2double Inlet_SpeciesVal_Sum = 0.0;
        for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          checkScalarBounds(Inlet_SpeciesVal[iMarker][iSpecies], "MARKER_INLET_SPECIES individual", 0.0, 1.0);
          Inlet_SpeciesVal_Sum += Inlet_SpeciesVal[iMarker][iSpecies];
        }
        checkScalarBounds(Inlet_SpeciesVal_Sum, "MARKER_INLET_SPECIES sum", 0.0, 1.0);
      }
    }

  } // species transport checks

  /*--- Define some variables for flamelet model. ---*/
  if (Kind_Species_Model == SPECIES_MODEL::FLAMELET) {
    /*--- The controlling variables are progress variable, total enthalpy, and optionally mixture fraction ---*/
    if (n_control_vars != (nSpecies - n_user_scalars))
      SU2_MPI::Error("Number of initial species incompatible with number of controlling variables and user scalars.", CURRENT_FUNCTION);
    /*--- We can have additional user defined transported scalars ---*/
    n_scalars = n_control_vars + n_user_scalars;
  }

  if (Kind_Regime == ENUM_REGIME::COMPRESSIBLE && GetBounded_Scalar()) {
    SU2_MPI::Error("BOUNDED_SCALAR discretization can only be used for incompressible problems.", CURRENT_FUNCTION);
  }

}

void CConfig::SetMarkers(SU2_COMPONENT val_software) {

  unsigned short iMarker_All, iMarker_CfgFile, iMarker_Euler, iMarker_Custom,
  iMarker_FarField, iMarker_SymWall, iMarker_PerBound,
  iMarker_NearFieldBound, iMarker_Fluid_InterfaceBound,
  iMarker_Inlet, iMarker_Riemann, iMarker_Giles, iMarker_Outlet,
  iMarker_Smoluchowski_Maxwell,
  iMarker_Isothermal,iMarker_HeatFlux,iMarker_HeatTansfer,
  iMarker_EngineInflow, iMarker_EngineExhaust, iMarker_Damper,
  iMarker_Displacement, iMarker_Load, iMarker_Internal,
  iMarker_Monitoring, iMarker_Designing, iMarker_GeoEval, iMarker_Plotting, iMarker_Analyze,
  iMarker_DV, iMarker_Moving, iMarker_SobolevBC, iMarker_PyCustom, iMarker_Supersonic_Inlet, iMarker_Supersonic_Outlet,
  iMarker_Clamped, iMarker_ZoneInterface, iMarker_CHTInterface, iMarker_Load_Dir, iMarker_Disp_Dir,
  iMarker_Fluid_Load, iMarker_Deform_Mesh, iMarker_Deform_Mesh_Sym_Plane,
  iMarker_ActDiskInlet, iMarker_ActDiskOutlet,
  iMarker_Turbomachinery, iMarker_MixingPlaneInterface;

  int size = SINGLE_NODE;
  SU2_MPI::Comm_size(SU2_MPI::GetComm(), &size);

  /*--- Compute the total number of markers in the config file ---*/
  nMarker_CfgFile = nMarker_Euler + nMarker_FarField + nMarker_SymWall +
  nMarker_PerBound + nMarker_NearFieldBound + nMarker_Fluid_InterfaceBound +
  nMarker_CHTInterface + nMarker_Inlet + nMarker_Riemann + nMarker_Smoluchowski_Maxwell +
  nMarker_Giles + nMarker_Outlet + nMarker_Isothermal +
  nMarker_HeatFlux + nMarker_HeatTransfer +
  nMarker_EngineInflow + nMarker_EngineExhaust + nMarker_Internal +
  nMarker_Supersonic_Inlet + nMarker_Supersonic_Outlet + nMarker_Displacement + nMarker_Load +
  nMarker_Custom + nMarker_Damper + nMarker_Fluid_Load +
  nMarker_Clamped + nMarker_Load_Dir + nMarker_Disp_Dir +
  nMarker_ActDiskInlet + nMarker_ActDiskOutlet +
  nMarker_ActDiskBemInlet_CG + nMarker_ActDiskBemOutlet_CG +
  nMarker_ZoneInterface;

  /*--- Add the possible send/receive domains ---*/

  nMarker_Max = nMarker_CfgFile + OVERHEAD*size;

  /*--- Basic dimensionalization of the markers (worst scenario) ---*/

  nMarker_All = nMarker_Max;

  /*--- Allocate the memory (markers in each domain) ---*/

  Marker_All_TagBound       = new string[nMarker_All];    // Store the tag that correspond with each marker.
  Marker_All_SendRecv       = new short[nMarker_All] ();   // +#domain (send), -#domain (receive).
  Marker_All_KindBC         = new unsigned short[nMarker_All] (); // Store the kind of boundary condition.
  Marker_All_Monitoring     = new unsigned short[nMarker_All] (); // Store whether the boundary should be monitored.
  Marker_All_Designing      = new unsigned short[nMarker_All] (); // Store whether the boundary should be designed.
  Marker_All_Plotting       = new unsigned short[nMarker_All] (); // Store whether the boundary should be plotted.
  Marker_All_Analyze        = new unsigned short[nMarker_All] (); // Store whether the boundary should be plotted.
  Marker_All_ZoneInterface  = new unsigned short[nMarker_All] (); // Store whether the boundary is in the FSI interface.
  Marker_All_GeoEval        = new unsigned short[nMarker_All] (); // Store whether the boundary should be geometry evaluation.
  Marker_All_DV             = new unsigned short[nMarker_All] (); // Store whether the boundary should be affected by design variables.
  Marker_All_Moving         = new unsigned short[nMarker_All] (); // Store whether the boundary should be in motion.
  Marker_All_Deform_Mesh    = new unsigned short[nMarker_All] (); // Store whether the boundary is deformable.
  Marker_All_Deform_Mesh_Sym_Plane = new unsigned short[nMarker_All] (); //Store wheter the boundary will follow the deformation
  Marker_All_Fluid_Load     = new unsigned short[nMarker_All] (); // Store whether the boundary computes/applies fluid loads.
  Marker_All_PyCustom       = new unsigned short[nMarker_All] (); // Store whether the boundary is Python customizable.
  Marker_All_PerBound       = new short[nMarker_All] ();          // Store whether the boundary belongs to a periodic boundary.
  Marker_All_Turbomachinery       = new unsigned short[nMarker_All] (); // Store whether the boundary is in needed for Turbomachinery computations.
  Marker_All_TurbomachineryFlag   = new unsigned short[nMarker_All] (); // Store whether the boundary has a flag for Turbomachinery computations.
  Marker_All_MixingPlaneInterface = new unsigned short[nMarker_All] (); // Store whether the boundary has a in the MixingPlane interface.
  Marker_All_Giles                = new unsigned short[nMarker_All] (); // Store whether the boundary has is a Giles boundary.
  Marker_All_SobolevBC      = new unsigned short[nMarker_All] (); // Store wether the boundary should apply to the gradient smoothing.

  for (iMarker_All = 0; iMarker_All < nMarker_All; iMarker_All++) {
    Marker_All_TagBound[iMarker_All] = "SEND_RECEIVE";
  }

  /*--- Allocate the memory (markers in the config file) ---*/

  Marker_CfgFile_TagBound             = new string[nMarker_CfgFile];
  Marker_CfgFile_KindBC               = new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_Monitoring           = new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_Designing            = new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_Plotting             = new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_Analyze              = new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_GeoEval              = new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_ZoneInterface        = new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_DV                   = new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_Moving               = new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_Deform_Mesh          = new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_Deform_Mesh_Sym_Plane= new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_Fluid_Load           = new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_PerBound             = new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_Turbomachinery       = new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_TurbomachineryFlag   = new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_MixingPlaneInterface = new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_Giles                = new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_PyCustom             = new unsigned short[nMarker_CfgFile] ();
  Marker_CfgFile_SobolevBC            = new unsigned short[nMarker_CfgFile] ();

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = "SEND_RECEIVE";
  }

  /*--- Allocate memory to store surface information (Analyze BC) ---*/

  Surface_MassFlow = new su2double[nMarker_Analyze] ();
  Surface_Mach = new su2double[nMarker_Analyze] ();
  Surface_Temperature = new su2double[nMarker_Analyze] ();
  Surface_Pressure = new su2double[nMarker_Analyze] ();
  Surface_Density = new su2double[nMarker_Analyze] ();
  Surface_Enthalpy = new su2double[nMarker_Analyze] ();
  Surface_NormalVelocity = new su2double[nMarker_Analyze] ();
  Surface_Uniformity = new su2double[nMarker_Analyze] ();
  Surface_SecondaryStrength = new su2double[nMarker_Analyze] ();
  Surface_SecondOverUniform = new su2double[nMarker_Analyze] ();
  Surface_MomentumDistortion = new su2double[nMarker_Analyze] ();
  Surface_TotalTemperature = new su2double[nMarker_Analyze] ();
  Surface_TotalPressure = new su2double[nMarker_Analyze] ();
  Surface_PressureDrop = new su2double[nMarker_Analyze] ();
  Surface_Species_0 = new su2double[nMarker_Analyze] ();
  Surface_Species_Variance = new su2double[nMarker_Analyze] ();
  Surface_DC60 = new su2double[nMarker_Analyze] ();
  Surface_IDC = new su2double[nMarker_Analyze] ();
  Surface_IDC_Mach = new su2double[nMarker_Analyze] ();
  Surface_IDR = new su2double[nMarker_Analyze] ();

  /*--- Populate the marker information in the config file (all domains) ---*/

  iMarker_CfgFile = 0;
  for (iMarker_Euler = 0; iMarker_Euler < nMarker_Euler; iMarker_Euler++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Euler[iMarker_Euler];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = EULER_WALL;
    iMarker_CfgFile++;
  }

  for (iMarker_FarField = 0; iMarker_FarField < nMarker_FarField; iMarker_FarField++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_FarField[iMarker_FarField];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = FAR_FIELD;
    iMarker_CfgFile++;
  }

  for (iMarker_SymWall = 0; iMarker_SymWall < nMarker_SymWall; iMarker_SymWall++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_SymWall[iMarker_SymWall];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = SYMMETRY_PLANE;
    iMarker_CfgFile++;
  }

  for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_PerBound[iMarker_PerBound];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = PERIODIC_BOUNDARY;
    Marker_CfgFile_PerBound[iMarker_CfgFile] = iMarker_PerBound + 1;
    iMarker_CfgFile++;
  }

  ActDisk_DeltaPress = new su2double[nMarker_ActDiskInlet] ();
  ActDisk_DeltaTemp = new su2double[nMarker_ActDiskInlet] ();
  ActDisk_TotalPressRatio = new su2double[nMarker_ActDiskInlet] ();
  ActDisk_TotalTempRatio = new su2double[nMarker_ActDiskInlet] ();
  ActDisk_StaticPressRatio = new su2double[nMarker_ActDiskInlet] ();
  ActDisk_StaticTempRatio = new su2double[nMarker_ActDiskInlet] ();
  ActDisk_Power = new su2double[nMarker_ActDiskInlet] ();
  ActDisk_MassFlow = new su2double[nMarker_ActDiskInlet] ();
  ActDisk_Mach = new su2double[nMarker_ActDiskInlet] ();
  ActDisk_Force = new su2double[nMarker_ActDiskInlet] ();
  ActDisk_NetThrust = new su2double[nMarker_ActDiskInlet] ();
  ActDisk_BCThrust = new su2double[nMarker_ActDiskInlet] ();
  ActDisk_BCThrust_Old = new su2double[nMarker_ActDiskInlet] ();
  ActDisk_GrossThrust = new su2double[nMarker_ActDiskInlet] ();
  ActDisk_Area = new su2double[nMarker_ActDiskInlet] ();
  ActDisk_ReverseMassFlow = new su2double[nMarker_ActDiskInlet] ();

  ActDiskInlet_MassFlow = new su2double[nMarker_ActDiskInlet] ();
  ActDiskInlet_Temperature = new su2double[nMarker_ActDiskInlet] ();
  ActDiskInlet_TotalTemperature = new su2double[nMarker_ActDiskInlet] ();
  ActDiskInlet_Pressure = new su2double[nMarker_ActDiskInlet] ();
  ActDiskInlet_TotalPressure = new su2double[nMarker_ActDiskInlet] ();
  ActDiskInlet_RamDrag = new su2double[nMarker_ActDiskInlet] ();
  ActDiskInlet_Force = new su2double[nMarker_ActDiskInlet] ();
  ActDiskInlet_Power = new su2double[nMarker_ActDiskInlet] ();

  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_ActDiskInlet[iMarker_ActDiskInlet];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = ACTDISK_INLET;
    iMarker_CfgFile++;
  }

  ActDiskOutlet_MassFlow = new su2double[nMarker_ActDiskOutlet] ();
  ActDiskOutlet_Temperature = new su2double[nMarker_ActDiskOutlet] ();
  ActDiskOutlet_TotalTemperature = new su2double[nMarker_ActDiskOutlet] ();
  ActDiskOutlet_Pressure = new su2double[nMarker_ActDiskOutlet] ();
  ActDiskOutlet_TotalPressure = new su2double[nMarker_ActDiskOutlet] ();
  ActDiskOutlet_GrossThrust = new su2double[nMarker_ActDiskOutlet] ();
  ActDiskOutlet_Force = new su2double[nMarker_ActDiskOutlet] ();
  ActDiskOutlet_Power = new su2double[nMarker_ActDiskOutlet] ();

  ActDiskOutlet_Thrust_BEM = new su2double[nMarker_ActDiskOutlet]();
  ActDiskOutlet_Torque_BEM = new su2double[nMarker_ActDiskOutlet]();

  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_ActDiskOutlet[iMarker_ActDiskOutlet];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = ACTDISK_OUTLET;
    iMarker_CfgFile++;
  }

  Outlet_MassFlow = new su2double[nMarker_Outlet] ();
  Outlet_Density  = new su2double[nMarker_Outlet] ();
  Outlet_Area     = new su2double[nMarker_Outlet] ();

  for (iMarker_NearFieldBound = 0; iMarker_NearFieldBound < nMarker_NearFieldBound; iMarker_NearFieldBound++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_NearFieldBound[iMarker_NearFieldBound];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = NEARFIELD_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_Fluid_InterfaceBound = 0; iMarker_Fluid_InterfaceBound < nMarker_Fluid_InterfaceBound; iMarker_Fluid_InterfaceBound++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Fluid_InterfaceBound[iMarker_Fluid_InterfaceBound];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = BC_TYPE::FLUID_INTERFACE;
    iMarker_CfgFile++;
  }

  for (iMarker_CHTInterface = 0; iMarker_CHTInterface < nMarker_CHTInterface; iMarker_CHTInterface++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_CHTInterface[iMarker_CHTInterface];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = CHT_WALL_INTERFACE;
    iMarker_CfgFile++;
  }

  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Inlet[iMarker_Inlet];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = INLET_FLOW;
    iMarker_CfgFile++;
  }

  for (iMarker_Riemann = 0; iMarker_Riemann < nMarker_Riemann; iMarker_Riemann++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Riemann[iMarker_Riemann];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = RIEMANN_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_Giles = 0; iMarker_Giles < nMarker_Giles; iMarker_Giles++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Giles[iMarker_Giles];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = GILES_BOUNDARY;
    iMarker_CfgFile++;
  }

  Engine_Power       = new su2double[nMarker_EngineInflow] ();
  Engine_Mach        = new su2double[nMarker_EngineInflow] ();
  Engine_Force       = new su2double[nMarker_EngineInflow] ();
  Engine_NetThrust   = new su2double[nMarker_EngineInflow] ();
  Engine_GrossThrust = new su2double[nMarker_EngineInflow] ();
  Engine_Area        = new su2double[nMarker_EngineInflow] ();

  Inflow_Mach = new su2double[nMarker_EngineInflow] ();
  Inflow_Pressure = new su2double[nMarker_EngineInflow] ();
  Inflow_MassFlow = new su2double[nMarker_EngineInflow] ();
  Inflow_ReverseMassFlow = new su2double[nMarker_EngineInflow] ();
  Inflow_TotalPressure = new su2double[nMarker_EngineInflow] ();
  Inflow_Temperature = new su2double[nMarker_EngineInflow] ();
  Inflow_TotalTemperature = new su2double[nMarker_EngineInflow] ();
  Inflow_RamDrag = new su2double[nMarker_EngineInflow] ();
  Inflow_Force = new su2double[nMarker_EngineInflow] ();
  Inflow_Power = new su2double[nMarker_EngineInflow] ();

  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_EngineInflow[iMarker_EngineInflow];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = ENGINE_INFLOW;
    iMarker_CfgFile++;
  }

  Exhaust_Pressure = new su2double[nMarker_EngineExhaust] ();
  Exhaust_Temperature = new su2double[nMarker_EngineExhaust] ();
  Exhaust_MassFlow = new su2double[nMarker_EngineExhaust] ();
  Exhaust_TotalPressure = new su2double[nMarker_EngineExhaust] ();
  Exhaust_TotalTemperature = new su2double[nMarker_EngineExhaust] ();
  Exhaust_GrossThrust = new su2double[nMarker_EngineExhaust] ();
  Exhaust_Force = new su2double[nMarker_EngineExhaust] ();
  Exhaust_Power = new su2double[nMarker_EngineExhaust] ();

  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_EngineExhaust[iMarker_EngineExhaust];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = ENGINE_EXHAUST;
    iMarker_CfgFile++;
  }

  for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = SUPERSONIC_INLET;
    iMarker_CfgFile++;
  }

  for (iMarker_Supersonic_Outlet = 0; iMarker_Supersonic_Outlet < nMarker_Supersonic_Outlet; iMarker_Supersonic_Outlet++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Supersonic_Outlet[iMarker_Supersonic_Outlet];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = SUPERSONIC_OUTLET;
    iMarker_CfgFile++;
  }

  for (iMarker_Internal = 0; iMarker_Internal < nMarker_Internal; iMarker_Internal++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Internal[iMarker_Internal];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = INTERNAL_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_Custom = 0; iMarker_Custom < nMarker_Custom; iMarker_Custom++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Custom[iMarker_Custom];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = CUSTOM_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Outlet[iMarker_Outlet];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = OUTLET_FLOW;
    iMarker_CfgFile++;
  }

  for (iMarker_Isothermal = 0; iMarker_Isothermal < nMarker_Isothermal; iMarker_Isothermal++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Isothermal[iMarker_Isothermal];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = ISOTHERMAL;
    iMarker_CfgFile++;
  }

  for (iMarker_Smoluchowski_Maxwell = 0; iMarker_Smoluchowski_Maxwell < nMarker_Smoluchowski_Maxwell; iMarker_Smoluchowski_Maxwell++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Smoluchowski_Maxwell[iMarker_Smoluchowski_Maxwell];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = SMOLUCHOWSKI_MAXWELL;
    iMarker_CfgFile++;
  }

  for (iMarker_HeatFlux = 0; iMarker_HeatFlux < nMarker_HeatFlux; iMarker_HeatFlux++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_HeatFlux[iMarker_HeatFlux];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = HEAT_FLUX;
    iMarker_CfgFile++;
  }

  for (iMarker_HeatTansfer = 0; iMarker_HeatTansfer < nMarker_HeatTransfer; iMarker_HeatTansfer++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_HeatTransfer[iMarker_HeatTansfer];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = HEAT_TRANSFER;
    iMarker_CfgFile++;
  }

  for (iMarker_Clamped = 0; iMarker_Clamped < nMarker_Clamped; iMarker_Clamped++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Clamped[iMarker_Clamped];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = CLAMPED_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_Displacement = 0; iMarker_Displacement < nMarker_Displacement; iMarker_Displacement++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Displacement[iMarker_Displacement];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = DISPLACEMENT_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_Load = 0; iMarker_Load < nMarker_Load; iMarker_Load++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Load[iMarker_Load];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = LOAD_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_Damper = 0; iMarker_Damper < nMarker_Damper; iMarker_Damper++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Damper[iMarker_Damper];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = DAMPER_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_Load_Dir = 0; iMarker_Load_Dir < nMarker_Load_Dir; iMarker_Load_Dir++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Load_Dir[iMarker_Load_Dir];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = LOAD_DIR_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_Disp_Dir = 0; iMarker_Disp_Dir < nMarker_Disp_Dir; iMarker_Disp_Dir++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Disp_Dir[iMarker_Disp_Dir];
    Marker_CfgFile_KindBC[iMarker_CfgFile] = DISP_DIR_BOUNDARY;
    iMarker_CfgFile++;
  }

  for (iMarker_Fluid_Load = 0; iMarker_Fluid_Load < nMarker_Fluid_Load; iMarker_Fluid_Load++) {
    Marker_CfgFile_TagBound[iMarker_CfgFile] = Marker_Fluid_Load[iMarker_Fluid_Load];
    iMarker_CfgFile++;
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_Monitoring[iMarker_CfgFile] = NO;
    for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarker_Monitoring; iMarker_Monitoring++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_Monitoring[iMarker_Monitoring])
        Marker_CfgFile_Monitoring[iMarker_CfgFile] = YES;
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_GeoEval[iMarker_CfgFile] = NO;
    for (iMarker_GeoEval = 0; iMarker_GeoEval < nMarker_GeoEval; iMarker_GeoEval++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_GeoEval[iMarker_GeoEval])
        Marker_CfgFile_GeoEval[iMarker_CfgFile] = YES;
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_Designing[iMarker_CfgFile] = NO;
    for (iMarker_Designing = 0; iMarker_Designing < nMarker_Designing; iMarker_Designing++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_Designing[iMarker_Designing])
        Marker_CfgFile_Designing[iMarker_CfgFile] = YES;
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_Plotting[iMarker_CfgFile] = NO;
    for (iMarker_Plotting = 0; iMarker_Plotting < nMarker_Plotting; iMarker_Plotting++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_Plotting[iMarker_Plotting])
        Marker_CfgFile_Plotting[iMarker_CfgFile] = YES;
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_Analyze[iMarker_CfgFile] = NO;
    for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_Analyze[iMarker_Analyze])
        Marker_CfgFile_Analyze[iMarker_CfgFile] = YES;
  }

  /*--- Identification of multi-physics interface markers ---*/

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_ZoneInterface[iMarker_CfgFile] = NO;
    for (iMarker_ZoneInterface = 0; iMarker_ZoneInterface < nMarker_ZoneInterface; iMarker_ZoneInterface++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_ZoneInterface[iMarker_ZoneInterface])
        Marker_CfgFile_ZoneInterface[iMarker_CfgFile] = YES;
  }

  /*--- Identification of Turbomachinery markers and flag them---*/

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    unsigned short indexMarker=0;
    Marker_CfgFile_Turbomachinery[iMarker_CfgFile] = NO;
    Marker_CfgFile_TurbomachineryFlag[iMarker_CfgFile] = NO;
    for (iMarker_Turbomachinery = 0; iMarker_Turbomachinery < nMarker_Turbomachinery; iMarker_Turbomachinery++){
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_TurboBoundIn[iMarker_Turbomachinery]){
        indexMarker=(iMarker_Turbomachinery+1);
        Marker_CfgFile_Turbomachinery[iMarker_CfgFile] = indexMarker;
        Marker_CfgFile_TurbomachineryFlag[iMarker_CfgFile] = INFLOW;
      }
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_TurboBoundOut[iMarker_Turbomachinery]){
        indexMarker=(iMarker_Turbomachinery+1);
        Marker_CfgFile_Turbomachinery[iMarker_CfgFile] = indexMarker;
        Marker_CfgFile_TurbomachineryFlag[iMarker_CfgFile] = OUTFLOW;
      }
    }
  }

  /*--- Idenftification fo Giles Markers ---*/
  // This is seperate from MP and Turbomachinery Markers as all mixing plane markers are Giles,
  // but not all Giles markers are mixing plane
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_Giles[iMarker_CfgFile] = NO;
    for (iMarker_Giles = 0; iMarker_Giles < nMarker_Giles; iMarker_Giles++) {
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_Giles[iMarker_Giles])
        Marker_CfgFile_Giles[iMarker_CfgFile] = YES;
    }
  }

  /*--- Identification of MixingPlane interface markers ---*/

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    unsigned short indexMarker=0;
    Marker_CfgFile_MixingPlaneInterface[iMarker_CfgFile] = NO;
    for (iMarker_MixingPlaneInterface = 0; iMarker_MixingPlaneInterface < nMarker_MixingPlaneInterface; iMarker_MixingPlaneInterface++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_MixingPlaneInterface[iMarker_MixingPlaneInterface])
        indexMarker=(int)(iMarker_MixingPlaneInterface/2+1);
    Marker_CfgFile_MixingPlaneInterface[iMarker_CfgFile] = indexMarker;
  }

  /*--- Once we have identified the MixingPlane and Turbomachinery markers
   *    we next need to determine what type of interface between zones is
   *    used. It is convenient to do this here as it tidies up the interface
   *    preproccesing in CDriver ---*/
  if (nMarker_Turbomachinery != 0) {
    nTurboInterfaces = (nMarker_Turbomachinery -  1)*2; //Two markers per zone minus inlet & outlet
    Kind_TurboInterface.resize(nTurboInterfaces);
    /*--- Loop over all markers ---*/
    for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
      /*--- Identify mixing plane markers ---*/
      if (Marker_MixingPlaneInterface != nullptr){ // Necessary in cases where no mixing plane interfaces are defined
        if (Marker_CfgFile_MixingPlaneInterface[iMarker_CfgFile] != 0) { //Is a mixing plane
          /*--- Find which list position this marker is in turbomachinery markers ---*/
          const auto* target = std::find(Marker_Turbomachinery, &Marker_Turbomachinery[nMarker_Turbomachinery*2-1], Marker_CfgFile_TagBound[iMarker_CfgFile]);
          const auto target_index = target - Marker_Turbomachinery;
          /*--- Assert that we find the marker within the turbomachienry markers ---*/
          assert(target != &Marker_Turbomachinery[nMarker_Turbomachinery*2-1]);
          /*--- Assign the correct interface ---*/
          Kind_TurboInterface[target_index-1] = TURBO_INTERFACE_KIND::MIXING_PLANE; // Need to subtract 1 from index as to not consider the inlet an interface
        }
      }
      if (Marker_Fluid_InterfaceBound != nullptr){ // No fluid interfaces are defined in the config file (nullptr if no interfaces defined)
        if (Marker_CfgFile_KindBC[iMarker_CfgFile] == BC_TYPE::FLUID_INTERFACE) { // iMarker_CfgFile is a fluid interface
          const auto* target = std::find(Marker_Turbomachinery, &Marker_Turbomachinery[nMarker_Turbomachinery*2-1], Marker_CfgFile_TagBound[iMarker_CfgFile]);
          const auto target_index = target - Marker_Turbomachinery;
          Kind_TurboInterface[target_index-1] = TURBO_INTERFACE_KIND::FROZEN_ROTOR;
        }
      }
    }
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_DV[iMarker_CfgFile] = NO;
    for (iMarker_DV = 0; iMarker_DV < nMarker_DV; iMarker_DV++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_DV[iMarker_DV])
        Marker_CfgFile_DV[iMarker_CfgFile] = YES;
  }

  /*--- Add an extra check for DV_MARKER to make sure that any given marker
   *    name is recognized as an existing boundary in the problem. ---*/

  for (iMarker_DV = 0; iMarker_DV < nMarker_DV; iMarker_DV++) {
    bool found = false;
    for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_DV[iMarker_DV]) {
        found = true;
        break;
      }
    }

    if(!found) {
      if (nZone==1)
        SU2_MPI::Error("DV_MARKER contains marker names that do not exist in the lists of BCs in the config file.", CURRENT_FUNCTION);
      // In case of multiple zones, the markers might appear only in zonal config and not in the Master.
      // A loop over all zones would need to be included which is not straight forward as this can only be
      // checked once all zonal configs are read.
      else if (rank == MASTER_NODE)
        cout << "Warning: DV_MARKER contains marker names that do not exist in the lists of BCs of the master config file.\n"
                "Make sure the marker names exist in the zonal config files" << endl;
    }
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_Moving[iMarker_CfgFile] = NO;
    for (iMarker_Moving = 0; iMarker_Moving < nMarker_Moving; iMarker_Moving++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_Moving[iMarker_Moving])
        Marker_CfgFile_Moving[iMarker_CfgFile] = YES;
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_Deform_Mesh[iMarker_CfgFile] = NO;
    for (iMarker_Deform_Mesh = 0; iMarker_Deform_Mesh < nMarker_Deform_Mesh; iMarker_Deform_Mesh++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_Deform_Mesh[iMarker_Deform_Mesh])
        Marker_CfgFile_Deform_Mesh[iMarker_CfgFile] = YES;
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_Deform_Mesh_Sym_Plane[iMarker_CfgFile] = NO;
    for (iMarker_Deform_Mesh_Sym_Plane = 0; iMarker_Deform_Mesh_Sym_Plane < nMarker_Deform_Mesh_Sym_Plane; iMarker_Deform_Mesh_Sym_Plane++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_Deform_Mesh_Sym_Plane[iMarker_Deform_Mesh_Sym_Plane])
        Marker_CfgFile_Deform_Mesh_Sym_Plane[iMarker_CfgFile] = YES;
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_Fluid_Load[iMarker_CfgFile] = NO;
    for (iMarker_Fluid_Load = 0; iMarker_Fluid_Load < nMarker_Fluid_Load; iMarker_Fluid_Load++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_Fluid_Load[iMarker_Fluid_Load])
        Marker_CfgFile_Fluid_Load[iMarker_CfgFile] = YES;
  }

  for (iMarker_CfgFile=0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_PyCustom[iMarker_CfgFile] = NO;
    for(iMarker_PyCustom=0; iMarker_PyCustom < nMarker_PyCustom; iMarker_PyCustom++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_PyCustom[iMarker_PyCustom])
        Marker_CfgFile_PyCustom[iMarker_CfgFile] = YES;
  }

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    Marker_CfgFile_SobolevBC[iMarker_CfgFile] = NO;
    for (iMarker_SobolevBC = 0; iMarker_SobolevBC < nMarker_SobolevBC; iMarker_SobolevBC++)
      if (Marker_CfgFile_TagBound[iMarker_CfgFile] == Marker_SobolevBC[iMarker_SobolevBC])
        Marker_CfgFile_SobolevBC[iMarker_CfgFile] = YES;
  }

}

void CConfig::SetOutput(SU2_COMPONENT val_software, unsigned short val_izone) {

  unsigned short iMarker_Euler, iMarker_Custom, iMarker_FarField,
  iMarker_SymWall, iMarker_PerBound, iMarker_NearFieldBound,
  iMarker_Fluid_InterfaceBound, iMarker_Inlet, iMarker_Riemann,
  iMarker_Deform_Mesh, iMarker_Deform_Mesh_Sym_Plane, iMarker_Fluid_Load,
  iMarker_Smoluchowski_Maxwell, iWall_Catalytic,
  iMarker_Giles, iMarker_Outlet, iMarker_Isothermal, iMarker_HeatFlux, iMarker_HeatTransfer,
  iMarker_EngineInflow, iMarker_EngineExhaust, iMarker_Displacement, iMarker_Damper,
  iMarker_Load, iMarker_Internal, iMarker_Monitoring,
  iMarker_Designing, iMarker_GeoEval, iMarker_Plotting, iMarker_Analyze, iMarker_DV, iDV_Value,
  iMarker_ZoneInterface, iMarker_PyCustom, iMarker_Load_Dir, iMarker_Disp_Dir, iMarker_Clamped,
  iMarker_Moving, iMarker_Supersonic_Inlet, iMarker_Supersonic_Outlet, iMarker_ActDiskInlet,
  iMarker_Emissivity, iMarker_StrongBC,
  iMarker_ActDiskOutlet, iMarker_MixingPlaneInterface,
  iMarker_SobolevBC;

  bool fea = ((Kind_Solver == MAIN_SOLVER::FEM_ELASTICITY) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_FEM));

  cout << endl <<"----------------- Physical Case Definition ( Zone "  << iZone << " ) -------------------" << endl;
  if (val_software == SU2_COMPONENT::SU2_CFD) {
    if (FSI_Problem)
     cout << "Fluid-Structure Interaction." << endl;

    if (DiscreteAdjoint) {
     cout <<"Discrete Adjoint equations using Algorithmic Differentiation\n";
     cout <<"based on the physical case: ";
    }
    switch (Kind_Solver) {
      case MAIN_SOLVER::EULER:     case MAIN_SOLVER::DISC_ADJ_EULER:
      case MAIN_SOLVER::INC_EULER: case MAIN_SOLVER::DISC_ADJ_INC_EULER:
      case MAIN_SOLVER::FEM_EULER: case MAIN_SOLVER::DISC_ADJ_FEM_EULER:
        if (Kind_Regime == ENUM_REGIME::COMPRESSIBLE) cout << "Compressible Euler equations." << endl;
        if (Kind_Regime == ENUM_REGIME::INCOMPRESSIBLE) cout << "Incompressible Euler equations." << endl;
        break;
      case MAIN_SOLVER::NAVIER_STOKES:     case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES:
      case MAIN_SOLVER::INC_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES:
      case MAIN_SOLVER::FEM_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_FEM_NS:
        if (Kind_Regime == ENUM_REGIME::COMPRESSIBLE) cout << "Compressible Laminar Navier-Stokes' equations." << endl;
        if (Kind_Regime == ENUM_REGIME::INCOMPRESSIBLE) cout << "Incompressible Laminar Navier-Stokes' equations." << endl;
        break;
      case MAIN_SOLVER::RANS:     case MAIN_SOLVER::DISC_ADJ_RANS:
      case MAIN_SOLVER::INC_RANS: case MAIN_SOLVER::DISC_ADJ_INC_RANS:
      case MAIN_SOLVER::FEM_RANS: case MAIN_SOLVER::DISC_ADJ_FEM_RANS:
        if (Kind_Regime == ENUM_REGIME::COMPRESSIBLE) cout << "Compressible RANS equations." << endl;
        if (Kind_Regime == ENUM_REGIME::INCOMPRESSIBLE) cout << "Incompressible RANS equations." << endl;
        cout << "Turbulence model: ";
        switch (Kind_Turb_Model) {
          case TURB_MODEL::NONE: break;
          case TURB_MODEL::SA:
            switch (saParsedOptions.version) {
              case SA_OPTIONS::NEG:
                cout << "Negative-";
                break;
              case SA_OPTIONS::EDW:
                cout << "Edwards-";
                break;
              default:
                break;
            }
            cout << "Spalart-Allmaras";

            if (!saParsedOptions.ft2) cout << "-noft2";
            if (saParsedOptions.rot) cout << "-R";
            if (saParsedOptions.comp) cout << "-comp";
            if (saParsedOptions.qcr2000) cout << "-QCR2000";
            if (saParsedOptions.bc) cout << "-BCM";
            cout << endl;
            break;
          case TURB_MODEL::SST:
            cout << "Menter's k-omega SST";
            if (sstParsedOptions.version == SST_OPTIONS::V1994) cout << "-1994";
            else cout << "-2003";
            if (sstParsedOptions.modified) cout << "m";
            if (sstParsedOptions.sust) cout << " with sustaining terms,";

            switch (sstParsedOptions.production) {
              case SST_OPTIONS::KL:
                cout << " with Kato-Launder production";
                break;
              case SST_OPTIONS::V:
                cout << " with Vorticity production";
                break;
              case SST_OPTIONS::UQ:
                cout << "\nperturbing the Reynold's Stress Matrix towards " << eig_val_comp << " component turbulence";
                if (uq_permute) cout << " (permuting eigenvectors)";
                break;
              case SST_OPTIONS::COMP_Wilcox:
                cout << " with compressibility correction of Wilcox";
                break;
              case SST_OPTIONS::COMP_Sarkar:
                cout << " with compressibility correction of Sarkar";
                break;
              default:
                cout << " with no production modification";
                break;
            }

            if (sstParsedOptions.dll){
              cout << "\nusing non dimensional lower limits relative to infinity values clipping by Coefficients:" ;
              cout << " C_w= " << OmegaFactor_LowerLimit << " and C_k= " <<KFactor_LowerLimit ;
            }
            else cout << "\nusing default hard coded lower limit clipping";

            cout << "." << endl;
            break;
        }
        switch (Kind_Trans_Model) {
          case TURB_TRANS_MODEL::NONE:  break;
          case TURB_TRANS_MODEL::LM: {
            cout << "Transition model: Langtry and Menter's 4 equation model";
            if (lmParsedOptions.LM2015) {
              cout << " w/ cross-flow corrections (2015)" << endl;
            } else {
              cout << " (2009)" << endl;
            }
            break;
          }
          case TURB_TRANS_MODEL::AFT:{
            cout << "Transition model: Amplification Factor Transport model";
            break;
          }
        }
        if (Kind_Trans_Model == TURB_TRANS_MODEL::LM) {

          cout << "Correlation Functions: ";
          switch (lmParsedOptions.Correlation) {
            case TURB_TRANS_CORRELATION::MALAN: cout << "Malan et al. (2009)" << endl;  break;
            case TURB_TRANS_CORRELATION::SULUKSNA: cout << "Suluksna et al. (2009)" << endl;  break;
            case TURB_TRANS_CORRELATION::KRAUSE: cout << "Krause et al. (2008)" << endl;  break;
            case TURB_TRANS_CORRELATION::KRAUSE_HYPER: cout << "Krause et al. (2008, paper)" << endl;  break;
            case TURB_TRANS_CORRELATION::MEDIDA_BAEDER: cout << "Medida and Baeder (2011)" << endl;  break;
            case TURB_TRANS_CORRELATION::MEDIDA: cout << "Medida PhD (2014)" << endl;  break;
            case TURB_TRANS_CORRELATION::MENTER_LANGTRY: cout << "Menter and Langtry (2009)" << endl;  break;
            case TURB_TRANS_CORRELATION::DEFAULT:
              switch (Kind_Turb_Model) {
                case TURB_MODEL::SA: cout << "Malan et al. (2009)" << endl;  break;
                case TURB_MODEL::SST: cout << "Menter and Langtry (2009)" << endl;  break;
                case TURB_MODEL::NONE: SU2_MPI::Error("No turbulence model has been selected but LM transition model is active.", CURRENT_FUNCTION); break;
              }
              break;
          }
        }
        if (Kind_Trans_Model == TURB_TRANS_MODEL::AFT) {

          switch (aftParsedOptions.Correlation) {
            case AFT_CORRELATION::AFT2017b:
              switch (Kind_Turb_Model) {
                case TURB_MODEL::NONE: SU2_MPI::Error("No turbulence model has been selected but AFT transition model is active.", CURRENT_FUNCTION); break;
                case TURB_MODEL::SST: SU2_MPI::Error("k-w SST turbulence model has been selected but AFT transition model is active.", CURRENT_FUNCTION); break;
              }
              cout << "-2017b" << endl;  break;
              if(!saParsedOptions.ft2) SU2_MPI::Error("ft2 option of SA model has been not selected.", CURRENT_FUNCTION);
            case AFT_CORRELATION::AFT2019b:
              switch (Kind_Turb_Model) {
                case TURB_MODEL::NONE: SU2_MPI::Error("No turbulence model has been selected but AFT transition model is active.", CURRENT_FUNCTION); break;
                case TURB_MODEL::SST: SU2_MPI::Error("k-w SST turbulence model has been selected but AFT transition model is active.", CURRENT_FUNCTION); break;
              }
              cout << "-2019b" << endl;  break;
              if(!saParsedOptions.ft2) SU2_MPI::Error("ft2 option of SA model has been not selected.", CURRENT_FUNCTION);
            case AFT_CORRELATION::NONE: SU2_MPI::Error("NONE has been selected.", CURRENT_FUNCTION); break;
          }
        }
        cout << "Hybrid RANS/LES: ";
        switch (Kind_HybridRANSLES) {
          case NO_HYBRIDRANSLES: cout << "No Hybrid RANS/LES" << endl; break;
          case SA_DES:   cout << "Detached Eddy Simulation (DES97) " << endl; break;
          case SA_DDES:  cout << "Delayed Detached Eddy Simulation (DDES) with Standard SGS" << endl; break;
          case SA_ZDES:  cout << "Delayed Detached Eddy Simulation (DDES) with Vorticity-based SGS" << endl; break;
          case SA_EDDES: cout << "Delayed Detached Eddy Simulation (DDES) with Shear-layer Adapted SGS" << endl; break;
        }
        break;
      case MAIN_SOLVER::NEMO_EULER:
        if (Kind_Regime == ENUM_REGIME::COMPRESSIBLE) cout << "Compressible two-temperature thermochemical non-equilibrium Euler equations." << endl;
        if (Kind_FluidModel == SU2_NONEQ){
          if ((GasModel != "N2") && (GasModel != "AIR-5") && (GasModel != "AIR-7") && (GasModel != "ARGON"))
            SU2_MPI::Error("The GAS_MODEL given is unavailable using CSU2TCLIB. Choose one of the options: N2, AIR-5, AIR-7, or ARGON.", CURRENT_FUNCTION);
        }
        break;
      case MAIN_SOLVER::NEMO_NAVIER_STOKES:
        if (Kind_Regime == ENUM_REGIME::COMPRESSIBLE) cout << "Compressible two-temperature thermochemical non-equilibrium Navier-Stokes equations." << endl;
        if (Kind_FluidModel == SU2_NONEQ){
          if ((GasModel != "N2") && (GasModel != "AIR-5") && (GasModel != "AIR-7") && (GasModel != "ARGON"))
          SU2_MPI::Error("The GAS_MODEL given is unavailable using CSU2TCLIB. Choose one of the options: N2, AIR-5, AIR-7, or ARGON.", CURRENT_FUNCTION);
        }
        break;
      case MAIN_SOLVER::FEM_LES:
        if (Kind_Regime == ENUM_REGIME::COMPRESSIBLE)   cout << "Compressible LES equations." << endl;
        if (Kind_Regime == ENUM_REGIME::INCOMPRESSIBLE) cout << "Incompressible LES equations." << endl;
        cout << "LES Subgrid Scale model: ";
        switch (Kind_SGS_Model) {
          case TURB_SGS_MODEL::IMPLICIT_LES: cout << "Implicit LES" << endl; break;
          case TURB_SGS_MODEL::SMAGORINSKY:  cout << "Smagorinsky " << endl; break;
          case TURB_SGS_MODEL::WALE:         cout << "WALE"         << endl; break;
          case TURB_SGS_MODEL::VREMAN:       cout << "VREMAN"         << endl; break;
          default:
            SU2_MPI::Error("Subgrid Scale model not specified.", CURRENT_FUNCTION);

        }
        break;
      case MAIN_SOLVER::FEM_ELASTICITY: case MAIN_SOLVER::DISC_ADJ_FEM:
        if (Kind_Struct_Solver == STRUCT_DEFORMATION::SMALL) cout << "Geometrically linear elasticity solver." << endl;
        if (Kind_Struct_Solver == STRUCT_DEFORMATION::LARGE) cout << "Geometrically non-linear elasticity solver." << endl;
        if (Kind_Material == STRUCT_MODEL::LINEAR_ELASTIC) cout << "Linear elastic material." << endl;
        if (Kind_Material == STRUCT_MODEL::NEO_HOOKEAN) {
          if (Kind_Material_Compress == STRUCT_COMPRESS::COMPRESSIBLE)
            cout << "Compressible Neo-Hookean material model." << endl;
        }
        break;
      case MAIN_SOLVER::ADJ_EULER: cout << "Continuous Euler adjoint equations." << endl; break;
      case MAIN_SOLVER::ADJ_NAVIER_STOKES:
        if (Frozen_Visc_Cont)
          cout << "Continuous Navier-Stokes adjoint equations with frozen (laminar) viscosity." << endl;
        else
          cout << "Continuous Navier-Stokes adjoint equations." << endl;
        break;
      case MAIN_SOLVER::ADJ_RANS:
        if (Frozen_Visc_Cont)
          cout << "Continuous RANS adjoint equations with frozen (laminar and eddy) viscosity." << endl;
        else
          cout << "Continuous RANS adjoint equations." << endl;
        break;
      case MAIN_SOLVER::HEAT_EQUATION: case MAIN_SOLVER::DISC_ADJ_HEAT:
        cout << "Heat solver" << endl;
        break;
      case MAIN_SOLVER::MULTIPHYSICS:
        cout << "Multiphysics solver" << endl;
        break;
      default:
        SU2_MPI::Error("No valid solver was chosen", CURRENT_FUNCTION);

    }

    if ((Kind_Regime == ENUM_REGIME::COMPRESSIBLE) && (Kind_Solver != MAIN_SOLVER::FEM_ELASTICITY)) {
      cout << "Mach number: " << Mach <<"."<< endl;
      cout << "Angle of attack (AoA): " << AoA <<" deg, and angle of sideslip (AoS): " << AoS <<" deg."<< endl;
      if ((Kind_Solver == MAIN_SOLVER::NAVIER_STOKES) || (Kind_Solver == MAIN_SOLVER::ADJ_NAVIER_STOKES) ||
          (Kind_Solver == MAIN_SOLVER::RANS) || (Kind_Solver == MAIN_SOLVER::ADJ_RANS) ||
          (Kind_Solver == MAIN_SOLVER::NEMO_NAVIER_STOKES))
        cout << "Reynolds number: " << Reynolds <<". Reference length "  << Length_Reynolds << "." << endl;
      if (Fixed_CL_Mode) {
        cout << "Fixed CL mode, target value: " << Target_CL << "." << endl;
      }
    }

    if (EquivArea) {
      cout <<"The equivalent area is going to be evaluated on the near-field."<< endl;
      cout <<"The lower integration limit is "<<ea_lim[0]<<", and the upper is "<<ea_lim[1]<<"."<< endl;
      cout <<"The near-field is situated at "<<ea_lim[2]<<"."<< endl;
    }

    if (GetGrid_Movement()) {
      cout << "Performing a dynamic mesh simulation: ";
      switch (Kind_GridMovement) {
        case NO_MOVEMENT:     cout << "no direct movement." << endl; break;
        case RIGID_MOTION:    cout << "rigid mesh motion." << endl; break;
        case ROTATING_FRAME:  cout << "rotating reference frame." << endl; break;
        case EXTERNAL:        cout << "externally prescribed motion." << endl; break;
      }
    }

    if (Restart) {
      if (Read_Binary_Restart) cout << "Reading and writing binary SU2 native restart files." << endl;
      else cout << "Reading and writing ASCII SU2 native restart files." << endl;
      if (!ContinuousAdjoint && Kind_Solver != MAIN_SOLVER::FEM_ELASTICITY) cout << "Read flow solution from: " << Solution_FileName << "." << endl;
      if (ContinuousAdjoint) cout << "Read adjoint solution from: " << Solution_AdjFileName << "." << endl;
    }
    else {
        if (fea) cout << "No restart solution, initialize from undeformed configuration." << endl;
        else cout << "No restart solution, use the values at infinity (freestream)." << endl;
    }

    if (ContinuousAdjoint)
      cout << "Read flow solution from: " << Solution_FileName << "." << endl;

    if (!fea){
      if (Kind_Regime == ENUM_REGIME::COMPRESSIBLE) {
        if (Ref_NonDim == DIMENSIONAL) { cout << "Dimensional simulation." << endl; }
        else if (Ref_NonDim == FREESTREAM_PRESS_EQ_ONE) { cout << "Non-Dimensional simulation (P=1.0, Rho=1.0, T=1.0 at the farfield)." << endl; }
        else if (Ref_NonDim == FREESTREAM_VEL_EQ_MACH) { cout << "Non-Dimensional simulation (V=Mach, Rho=1.0, T=1.0 at the farfield)." << endl; }
        else if (Ref_NonDim == FREESTREAM_VEL_EQ_ONE) { cout << "Non-Dimensional simulation (V=1.0, Rho=1.0, T=1.0 at the farfield)." << endl; }
    } else if (Kind_Regime == ENUM_REGIME::INCOMPRESSIBLE) {
        if (Ref_Inc_NonDim == DIMENSIONAL) { cout << "Dimensional simulation." << endl; }
        else if (Ref_Inc_NonDim == INITIAL_VALUES) { cout << "Non-Dimensional simulation using intialization values." << endl; }
        else if (Ref_Inc_NonDim == REFERENCE_VALUES) { cout << "Non-Dimensional simulation using user-specified reference values." << endl; }
      }

      if (RefArea == 0.0) cout << "The reference area will be computed using y(2D) or z(3D) projection." << endl;
      else { cout << "The reference area is " << RefArea;
        if (SystemMeasurements == US) cout << " ft^2." << endl; else cout << " m^2." << endl;
      }

      if (SemiSpan == 0.0) cout << "The semi-span will be computed using the max y(3D) value." << endl;
      else { cout << "The semi-span length area is " << SemiSpan;
        if (SystemMeasurements == US) cout << " ft." << endl; else cout << " m." << endl;
      }

      cout << "The reference length is " << RefLength;
      if (SystemMeasurements == US) cout << " ft." << endl; else cout << " m." << endl;

      if (nMarker_Monitoring != 0){
        if ((nRefOriginMoment_X > 1) || (nRefOriginMoment_Y > 1) || (nRefOriginMoment_Z > 1)) {
          cout << "Surface(s) where the force coefficients are evaluated and \n";
          cout << "their reference origin for moment computation: \n";

          for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarker_Monitoring; iMarker_Monitoring++) {
            cout << "   - " << Marker_Monitoring[iMarker_Monitoring] << " (" << RefOriginMoment_X[iMarker_Monitoring] <<", "<<RefOriginMoment_Y[iMarker_Monitoring] <<", "<< RefOriginMoment_Z[iMarker_Monitoring] << ")";
            if (iMarker_Monitoring < nMarker_Monitoring-1) cout << ".\n";
            else {
              if (SystemMeasurements == US) cout <<" ft."<< endl;
              else cout <<" m."<< endl;
            }

          }
        }
        else {
          cout << "Reference origin for moment evaluation is (" << RefOriginMoment_X[0] << ", " << RefOriginMoment_Y[0] << ", " << RefOriginMoment_Z[0] << ")." << endl;
          cout << "Surface(s) where the force coefficients are evaluated: ";
          for (iMarker_Monitoring = 0; iMarker_Monitoring < nMarker_Monitoring; iMarker_Monitoring++) {
            cout << Marker_Monitoring[iMarker_Monitoring];
            if (iMarker_Monitoring < nMarker_Monitoring-1) cout << ", ";
            else cout <<"."<< endl;
          }
          cout<< endl;
        }
      }
    }

    if (nMarker_Designing != 0) {
      cout << "Surface(s) where the objective function is evaluated: ";
      for (iMarker_Designing = 0; iMarker_Designing < nMarker_Designing; iMarker_Designing++) {
        cout << Marker_Designing[iMarker_Designing];
        if (iMarker_Designing < nMarker_Designing-1) cout << ", ";
        else cout <<".";
      }
      cout<< endl;
    }

    if (nMarker_Plotting != 0) {
      cout << "Surface(s) plotted in the output file: ";
      for (iMarker_Plotting = 0; iMarker_Plotting < nMarker_Plotting; iMarker_Plotting++) {
        cout << Marker_Plotting[iMarker_Plotting];
        if (iMarker_Plotting < nMarker_Plotting-1) cout << ", ";
        else cout <<".";
      }
      cout<< endl;
    }

    if (nMarker_Analyze != 0) {
      cout << "Surface(s) to be analyzed in detail: ";
      for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
        cout << Marker_Analyze[iMarker_Analyze];
        if (iMarker_Analyze < nMarker_Analyze-1) cout << ", ";
        else cout <<".";
      }
      cout<< endl;
    }

    if (nMarker_ZoneInterface != 0) {
      cout << "Surface(s) acting as an interface among zones: ";
      for (iMarker_ZoneInterface = 0; iMarker_ZoneInterface < nMarker_ZoneInterface; iMarker_ZoneInterface++) {
        cout << Marker_ZoneInterface[iMarker_ZoneInterface];
        if (iMarker_ZoneInterface < nMarker_ZoneInterface-1) cout << ", ";
        else cout <<".";
      }
      cout<<endl;
    }

    if(nMarker_PyCustom != 0) {
      cout << "Surface(s) that are customizable in Python: ";
      for(iMarker_PyCustom=0; iMarker_PyCustom < nMarker_PyCustom; iMarker_PyCustom++){
        cout << Marker_PyCustom[iMarker_PyCustom];
        if (iMarker_PyCustom < nMarker_PyCustom-1) cout << ", ";
        else cout << ".";
      }
      cout << endl;
    }

    if (nMarker_DV != 0) {
      cout << "Surface(s) affected by the design variables: ";
      for (iMarker_DV = 0; iMarker_DV < nMarker_DV; iMarker_DV++) {
        cout << Marker_DV[iMarker_DV];
        if (iMarker_DV < nMarker_DV-1) cout << ", ";
        else cout <<".";
      }
      cout<< endl;
    }

    if (nMarker_Moving != 0) {
      cout << "Surface(s) in motion: ";
      for (iMarker_Moving = 0; iMarker_Moving < nMarker_Moving; iMarker_Moving++) {
        cout << Marker_Moving[iMarker_Moving];
        if (iMarker_Moving < nMarker_Moving-1) cout << ", ";
        else cout <<".";
      }
      cout<< endl;
    }

  }

  if (val_software == SU2_COMPONENT::SU2_GEO) {
    if (nMarker_GeoEval != 0) {
      cout << "Surface(s) where the geometrical based functions is evaluated: ";
      for (iMarker_GeoEval = 0; iMarker_GeoEval < nMarker_GeoEval; iMarker_GeoEval++) {
        cout << Marker_GeoEval[iMarker_GeoEval];
        if (iMarker_GeoEval < nMarker_GeoEval-1) cout << ", ";
        else cout <<".";
      }
      cout<< endl;
    }
  }

  cout << "Input mesh file name: " << Mesh_FileName << endl;

  if (val_software == SU2_COMPONENT::SU2_DOT) {
    if (DiscreteAdjoint) {
      cout << "Input sensitivity file name: " << GetObjFunc_Extension(Solution_AdjFileName) << "." << endl;
    }else {
    cout << "Input sensitivity file name: " << SurfAdjCoeff_FileName << "." << endl;
  }
  }

  if (val_software == SU2_COMPONENT::SU2_DEF) {
    cout << endl <<"---------------- Grid deformation parameters ( Zone "  << iZone << " )  ----------------" << endl;
    cout << "Grid deformation using a linear elasticity method." << endl;

    if (Hold_GridFixed == YES) cout << "Hold some regions of the mesh fixed (hardcode implementation)." << endl;
  }

  if (val_software == SU2_COMPONENT::SU2_DOT) {
  cout << endl <<"-------------- Surface deformation parameters ( Zone "  << iZone << " ) ----------------" << endl;
  }

  if (((val_software == SU2_COMPONENT::SU2_DEF) || (val_software == SU2_COMPONENT::SU2_DOT)) && (Design_Variable[0] != NO_DEFORMATION)) {

    for (unsigned short iDV = 0; iDV < nDV; iDV++) {


      if ((Design_Variable[iDV] != NO_DEFORMATION) &&
          (Design_Variable[iDV] != FFD_SETTING) &&
          (Design_Variable[iDV] != SCALE_GRID) &&
          (Design_Variable[iDV] != TRANSLATE_GRID) &&
          (Design_Variable[iDV] != ROTATE_GRID) &&
          (Design_Variable[iDV] != SURFACE_FILE)) {

        if (iDV == 0)
          cout << "Design variables definition (markers <-> value <-> param):" << endl;

        switch (Design_Variable[iDV]) {
          case FFD_CONTROL_POINT_2D:  cout << "FFD 2D (control point) <-> "; break;
          case FFD_CAMBER_2D:         cout << "FFD 2D (camber) <-> "; break;
          case FFD_THICKNESS_2D:      cout << "FFD 2D (thickness) <-> "; break;
          case HICKS_HENNE:           cout << "Hicks Henne <-> " ; break;
          case SURFACE_BUMP:          cout << "Surface bump <-> " ; break;
          case ANGLE_OF_ATTACK:       cout << "Angle of attack <-> " ; break;
          case CST:                   cout << "Kulfan parameter number (CST) <-> " ; break;
          case TRANSLATION:           cout << "Translation design variable."; break;
          case SCALE:                 cout << "Scale design variable."; break;
          case NACA_4DIGITS:          cout << "NACA four digits <-> "; break;
          case PARABOLIC:             cout << "Parabolic <-> "; break;
          case AIRFOIL:               cout << "Airfoil <-> "; break;
          case ROTATION:              cout << "Rotation <-> "; break;
          case FFD_CONTROL_POINT:     cout << "FFD (control point) <-> "; break;
          case FFD_NACELLE:           cout << "FFD (nacelle) <-> "; break;
          case FFD_GULL:              cout << "FFD (gull) <-> "; break;
          case FFD_TWIST:             cout << "FFD (twist) <-> "; break;
          case FFD_ROTATION:          cout << "FFD (rotation) <-> "; break;
          case FFD_CONTROL_SURFACE:   cout << "FFD (control surface) <-> "; break;
          case FFD_CAMBER:            cout << "FFD (camber) <-> "; break;
          case FFD_THICKNESS:         cout << "FFD (thickness) -> "; break;
          case FFD_ANGLE_OF_ATTACK:   cout << "FFD (angle of attack) <-> "; break;
        }

        for (iMarker_DV = 0; iMarker_DV < nMarker_DV; iMarker_DV++) {
          cout << Marker_DV[iMarker_DV];
          if (iMarker_DV < nMarker_DV-1) cout << ", ";
          else cout << " <-> ";
        }

        for (iDV_Value = 0; iDV_Value < nDV_Value[iDV]; iDV_Value++) {
          cout << DV_Value[iDV][iDV_Value];
          if (iDV_Value != nDV_Value[iDV]-1) cout << ", ";
        }
        cout << " <-> ";

        if ((Design_Variable[iDV] == NO_DEFORMATION) ||
            (Design_Variable[iDV] == FFD_SETTING) ||
            (Design_Variable[iDV] == SCALE) ) nParamDV = 0;
        if (Design_Variable[iDV] == ANGLE_OF_ATTACK) nParamDV = 1;
        if ((Design_Variable[iDV] == FFD_CAMBER_2D) ||
            (Design_Variable[iDV] == FFD_THICKNESS_2D) ||
            (Design_Variable[iDV] == HICKS_HENNE) ||
            (Design_Variable[iDV] == PARABOLIC) ||
            (Design_Variable[iDV] == AIRFOIL) ||
            (Design_Variable[iDV] == FFD_GULL) ||
            (Design_Variable[iDV] == FFD_ANGLE_OF_ATTACK) ) nParamDV = 2;
        if ((Design_Variable[iDV] ==  TRANSLATION) ||
            (Design_Variable[iDV] ==  NACA_4DIGITS) ||
            (Design_Variable[iDV] ==  CST) ||
            (Design_Variable[iDV] ==  SURFACE_BUMP) ||
            (Design_Variable[iDV] ==  FFD_CAMBER) ||
            (Design_Variable[iDV] ==  FFD_THICKNESS) ) nParamDV = 3;
        if (Design_Variable[iDV] == FFD_CONTROL_POINT_2D) nParamDV = 5;
        if (Design_Variable[iDV] == ROTATION) nParamDV = 6;
        if ((Design_Variable[iDV] ==  FFD_CONTROL_POINT) ||
            (Design_Variable[iDV] ==  FFD_ROTATION) ||
            (Design_Variable[iDV] ==  FFD_CONTROL_SURFACE) ) nParamDV = 7;
        if (Design_Variable[iDV] == FFD_TWIST) nParamDV = 8;

        for (unsigned short iParamDV = 0; iParamDV < nParamDV; iParamDV++) {

          if (iParamDV == 0) cout << "( ";

          if ((iParamDV == 0) &&
              ((Design_Variable[iDV] == NO_DEFORMATION) ||
               (Design_Variable[iDV] == FFD_SETTING) ||
               (Design_Variable[iDV] == FFD_ANGLE_OF_ATTACK) ||
               (Design_Variable[iDV] == FFD_CONTROL_POINT_2D) ||
               (Design_Variable[iDV] == FFD_CAMBER_2D) ||
               (Design_Variable[iDV] == FFD_THICKNESS_2D) ||
               (Design_Variable[iDV] == FFD_CONTROL_POINT) ||
               (Design_Variable[iDV] == FFD_NACELLE) ||
               (Design_Variable[iDV] == FFD_GULL) ||
               (Design_Variable[iDV] == FFD_TWIST) ||
               (Design_Variable[iDV] == FFD_ROTATION) ||
               (Design_Variable[iDV] == FFD_CONTROL_SURFACE) ||
               (Design_Variable[iDV] == FFD_CAMBER) ||
               (Design_Variable[iDV] == FFD_THICKNESS))) cout << FFDTag[iDV];
          else cout << ParamDV[iDV][iParamDV];

          if (iParamDV < nParamDV-1) cout << ", ";
          else cout <<" )"<< endl;

        }

      }

      else if (Design_Variable[iDV] == NO_DEFORMATION) {
        cout << "No deformation of the numerical grid. Just output .su2 file." << endl;
      }

      else if (Design_Variable[iDV] == SCALE_GRID) {
        nParamDV = 0;
        cout << "Scaling of the volume grid by a constant factor." << endl;
      }

      else if (Design_Variable[iDV] == TRANSLATE_GRID) {
        nParamDV = 3;
        cout << "Rigid translation of the volume grid." << endl;
      }

      else if (Design_Variable[iDV] == ROTATE_GRID) {
        nParamDV = 6;
        cout << "Rigid rotation of the volume grid." << endl;
      }

      else if (Design_Variable[iDV] == FFD_SETTING) {

        cout << "Setting the FFD box structure." << endl;
        cout << "FFD boxes definition (FFD tag <-> degree <-> coord):" << endl;

        for (unsigned short iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++) {

          cout << TagFFDBox[iFFDBox] << " <-> ";

          for (unsigned short iDegreeFFD = 0; iDegreeFFD < 3; iDegreeFFD++) {
            if (iDegreeFFD == 0) cout << "( ";
            cout << DegreeFFDBox[iFFDBox][iDegreeFFD];
            if (iDegreeFFD < 2) cout << ", ";
            else cout <<" )";
          }

          cout << " <-> ";

          for (unsigned short iCoordFFD = 0; iCoordFFD < 24; iCoordFFD++) {
            if (iCoordFFD == 0) cout << "( ";
            cout << CoordFFDBox[iFFDBox][iCoordFFD];
            if (iCoordFFD < 23) cout << ", ";
            else cout <<" )"<< endl;
          }

        }

      }

      else cout << endl;

    }
  }

  if (((val_software == SU2_COMPONENT::SU2_CFD) && ( ContinuousAdjoint || DiscreteAdjoint)) || (val_software == SU2_COMPONENT::SU2_DOT)) {

    cout << endl <<"---------------- Design problem definition  ( Zone "  << iZone << " ) ------------------" << endl;
    if (nObj==1) {
      switch (Kind_ObjFunc[0]) {
        case DRAG_COEFFICIENT:           cout << "CD objective function";
          if (Fixed_CL_Mode) {           cout << " using fixed CL mode, dCD/dCL = " << dCD_dCL << "." << endl; }
          else {                         cout << "." << endl; }
          break;
        case LIFT_COEFFICIENT:           cout << "CL objective function." << endl; break;
        case MOMENT_X_COEFFICIENT:       cout << "CMx objective function" << endl;
          if (Fixed_CL_Mode) {           cout << " using fixed CL mode, dCMx/dCL = " << dCMx_dCL << "." << endl; }
          else {                         cout << "." << endl; }
          break;
        case MOMENT_Y_COEFFICIENT:       cout << "CMy objective function" << endl;
          if (Fixed_CL_Mode) {           cout << " using fixed CL mode, dCMy/dCL = " << dCMy_dCL << "." << endl; }
          else {                         cout << "." << endl; }
          break;
        case MOMENT_Z_COEFFICIENT:       cout << "CMz objective function" << endl;
          if (Fixed_CL_Mode) {           cout << " using fixed CL mode, dCMz/dCL = " << dCMz_dCL << "." << endl; }
          else {                         cout << "." << endl; }
          break;
        case INVERSE_DESIGN_PRESSURE:    cout << "Inverse design (Cp) objective function." << endl; break;
        case INVERSE_DESIGN_HEATFLUX:    cout << "Inverse design (Heat Flux) objective function." << endl; break;
        case SIDEFORCE_COEFFICIENT:      cout << "Side force objective function." << endl; break;
        case EFFICIENCY:                 cout << "CL/CD objective function." << endl; break;
        case EQUIVALENT_AREA:            cout << "Equivalent area objective function. CD weight: " << WeightCd <<"."<< endl;  break;
        case NEARFIELD_PRESSURE:         cout << "Nearfield pressure objective function. CD weight: " << WeightCd <<"."<< endl;  break;
        case FORCE_X_COEFFICIENT:        cout << "X-force objective function." << endl; break;
        case FORCE_Y_COEFFICIENT:        cout << "Y-force objective function." << endl; break;
        case FORCE_Z_COEFFICIENT:        cout << "Z-force objective function." << endl; break;
        case THRUST_COEFFICIENT:         cout << "Thrust objective function." << endl; break;
        case TORQUE_COEFFICIENT:         cout << "Torque efficiency objective function." << endl; break;
        case TOTAL_HEATFLUX:             cout << "Total heat flux objective function." << endl; break;
        case MAXIMUM_HEATFLUX:           cout << "Maximum heat flux objective function." << endl; break;
        case FIGURE_OF_MERIT:            cout << "Rotor Figure of Merit objective function." << endl; break;
        case BUFFET_SENSOR:              cout << "Buffet sensor objective function." << endl; break;
        case SURFACE_TOTAL_PRESSURE:     cout << "Average total pressure objective function." << endl; break;
        case SURFACE_STATIC_PRESSURE:    cout << "Average static pressure objective function." << endl; break;
        case SURFACE_STATIC_TEMPERATURE: cout << "Average static temperature objective function." << endl; break;
        case SURFACE_MASSFLOW:           cout << "Mass flow rate objective function." << endl; break;
        case SURFACE_MACH:               cout << "Mach number objective function." << endl; break;
        case CUSTOM_OBJFUNC:             cout << "Custom objective function." << endl; break;
        case REFERENCE_GEOMETRY:         cout << "Target geometry objective function." << endl; break;
        case REFERENCE_NODE:             cout << "Target node displacement objective function." << endl; break;
        case VOLUME_FRACTION:            cout << "Volume fraction objective function." << endl; break;
        case TOPOL_DISCRETENESS:         cout << "Topology discreteness objective function." << endl; break;
        case TOPOL_COMPLIANCE:           cout << "Topology compliance objective function." << endl; break;
        case STRESS_PENALTY:             cout << "Stress penalty objective function." << endl; break;
      }
    }
    else {
      cout << "Weighted sum objective function." << endl;
    }

  }

  if (val_software == SU2_COMPONENT::SU2_CFD) {

    auto PrintLimiterInfo = [&](const LIMITER kind_limiter) {
      cout << "Second order integration in space, with slope limiter.\n";
      switch (kind_limiter) {
        case LIMITER::NONE:
          cout << "No slope-limiting method. " << endl;
          break;
        case LIMITER::VENKATAKRISHNAN:
          cout << "Venkatakrishnan slope-limiting method, with constant: " << Venkat_LimiterCoeff << ".\n";
          cout << "The reference element size is: " << RefElemLength << ". " << endl;
          break;
        case LIMITER::NISHIKAWA_R3:
          cout << "Nishikawa's R3 slope-limiting method, with constant: " << Venkat_LimiterCoeff << ".\n";
          cout << "The reference element size is: " << RefElemLength << ". " << endl;
          break;
        case LIMITER::NISHIKAWA_R4:
          cout << "Nishikawa's R4 slope-limiting method, with constant: " << Venkat_LimiterCoeff << ".\n";
          cout << "The reference element size is: " << RefElemLength << ". " << endl;
          break;
        case LIMITER::NISHIKAWA_R5:
          cout << "Nishikawa's R5 slope-limiting method, with constant: " << Venkat_LimiterCoeff << ".\n";
          cout << "The reference element size is: " << RefElemLength << ". " << endl;
          break;
        case LIMITER::VENKATAKRISHNAN_WANG:
          cout << "Venkatakrishnan-Wang slope-limiting method, with constant: " << Venkat_LimiterCoeff << "." << endl;
          break;
        case LIMITER::BARTH_JESPERSEN:
          cout << "Barth-Jespersen slope-limiting method." << endl;
          break;
        case LIMITER::VAN_ALBADA_EDGE:
          cout << "Van Albada slope-limiting method implemented by edges." << endl;
          break;
        case LIMITER::SHARP_EDGES:
          cout << "Sharp edges slope-limiting method, with constant: " << Venkat_LimiterCoeff << ".\n";
          cout << "The reference element size is: " << RefElemLength << ".\n";
          cout << "The reference sharp edge distance is: " << AdjSharp_LimiterCoeff*RefElemLength*Venkat_LimiterCoeff << "." << endl;
          break;
        case LIMITER::WALL_DISTANCE:
          cout << "Wall distance slope-limiting method, with constant: " << Venkat_LimiterCoeff << ".\n";
          cout << "The reference element size is: " << RefElemLength << ".\n";
          cout << "The reference wall distance is: " << AdjSharp_LimiterCoeff*RefElemLength*Venkat_LimiterCoeff << "." << endl;
          break;
        default:
          SU2_MPI::Error("Unknown or invalid limiter type.", CURRENT_FUNCTION);
          break;
      }
    };

    cout << endl <<"--------------- Space Numerical Integration ( Zone "  << iZone << " ) ------------------" << endl;

    if ((Kind_Solver == MAIN_SOLVER::EULER)          || (Kind_Solver == MAIN_SOLVER::NAVIER_STOKES)          || (Kind_Solver == MAIN_SOLVER::RANS) ||
        (Kind_Solver == MAIN_SOLVER::INC_EULER)      || (Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES)      || (Kind_Solver == MAIN_SOLVER::INC_RANS) ||
        (Kind_Solver == MAIN_SOLVER::NEMO_EULER)     || (Kind_Solver == MAIN_SOLVER::NEMO_NAVIER_STOKES)     ||
        (Kind_Solver == MAIN_SOLVER::DISC_ADJ_EULER) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_RANS) ) {

      if (Kind_ConvNumScheme_Flow == SPACE_CENTERED) {
        if (Kind_Centered_Flow == CENTERED::LAX) {
          cout << "Lax-Friedrich scheme (1st order in space) for the flow inviscid terms.\n";
          cout << "Lax viscous coefficients (1st): " << Kappa_1st_Flow << ".\n";
          cout << "First order integration." << endl;
        }
        else {
          cout << "Jameson-Schmidt-Turkel scheme (2nd order in space) for the flow inviscid terms.\n";
          cout << "JST viscous coefficients (2nd & 4th): " << Kappa_2nd_Flow << ", " << Kappa_4th_Flow << ".\n";
          cout << "The method includes a grid stretching correction (p = 0.3)."<< endl;
        }
      }

      if (Kind_ConvNumScheme_Flow == SPACE_UPWIND) {
        if (Kind_Upwind_Flow == UPWIND::ROE)    cout << "Roe (with entropy fix = "<< EntropyFix_Coeff <<") solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == UPWIND::TURKEL) cout << "Roe-Turkel solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == UPWIND::AUSM)   cout << "AUSM solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == UPWIND::HLLC)   cout << "HLLC solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == UPWIND::SW)     cout << "Steger-Warming solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == UPWIND::MSW)    cout << "Modified Steger-Warming solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == UPWIND::L2ROE)  cout << "L2ROE Low Mach ROE solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == UPWIND::LMROE)  cout << "Rieper Low Mach ROE solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == UPWIND::SLAU)   cout << "Simple Low-Dissipation AUSM solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == UPWIND::SLAU2)  cout << "Simple Low-Dissipation AUSM 2 solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == UPWIND::FDS)    cout << "Flux difference splitting (FDS) upwind scheme for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == UPWIND::AUSMPLUSUP)  cout << "AUSM+-up solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == UPWIND::AUSMPLUSUP2) cout << "AUSM+-up2 solver for the flow inviscid terms."<< endl;
        if (Kind_Upwind_Flow == UPWIND::AUSMPLUSM)  cout << "AUSM+M solver for the flow inviscid terms."<< endl;

        if (Kind_Solver == MAIN_SOLVER::EULER         || Kind_Solver == MAIN_SOLVER::DISC_ADJ_EULER ||
            Kind_Solver == MAIN_SOLVER::NAVIER_STOKES || Kind_Solver == MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES ||
            Kind_Solver == MAIN_SOLVER::RANS          || Kind_Solver == MAIN_SOLVER::DISC_ADJ_RANS) {
          switch (Kind_RoeLowDiss) {
            case NO_ROELOWDISS: cout << "Standard Roe without low-dissipation function."<< endl; break;
            case NTS: cout << "Roe with NTS low-dissipation function."<< endl; break;
            case FD: cout << "Roe with DDES's FD low-dissipation function."<< endl; break;
            case NTS_DUCROS: cout << "Roe with NTS low-dissipation function + Ducros shock sensor."<< endl; break;
            case FD_DUCROS: cout << "Roe with DDES's FD low-dissipation function + Ducros shock sensor."<< endl; break;
          }
        }

        if (MUSCL_Flow) {
          PrintLimiterInfo(Kind_SlopeLimit_Flow);
        } else {
          cout << "First order integration in space." << endl;
        }

      }

    }

    if ((Kind_Solver == MAIN_SOLVER::RANS) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_RANS)) {
      if (Kind_ConvNumScheme_Turb == SPACE_UPWIND) {
        if (Kind_Upwind_Turb == UPWIND::SCALAR_UPWIND) cout << "Scalar upwind solver for the turbulence model." << endl;
        if (MUSCL_Turb) {
          PrintLimiterInfo(Kind_SlopeLimit_Turb);
        } else {
          cout << "First order integration in space." << endl;
        }
      }
    }

    if ((Kind_Solver == MAIN_SOLVER::ADJ_EULER) || (Kind_Solver == MAIN_SOLVER::ADJ_NAVIER_STOKES) || (Kind_Solver == MAIN_SOLVER::ADJ_RANS)) {

      if (Kind_ConvNumScheme_AdjFlow == SPACE_CENTERED) {
        if (Kind_Centered_AdjFlow == CENTERED::JST) {
          cout << "Jameson-Schmidt-Turkel scheme for the adjoint inviscid terms."<< endl;
          cout << "JST viscous coefficients (1st, 2nd, & 4th): " << Kappa_1st_AdjFlow
          << ", " << Kappa_2nd_AdjFlow << ", " << Kappa_4th_AdjFlow <<"."<< endl;
          cout << "The method includes a grid stretching correction (p = 0.3)."<< endl;
          cout << "Second order integration." << endl;
        }
        if (Kind_Centered_AdjFlow == CENTERED::LAX) {
          cout << "Lax-Friedrich scheme for the adjoint inviscid terms."<< endl;
          cout << "First order integration." << endl;
        }
      }

      if (Kind_ConvNumScheme_AdjFlow == SPACE_UPWIND) {
        if (Kind_Upwind_AdjFlow == UPWIND::ROE) cout << "Roe (with entropy fix = "<< EntropyFix_Coeff <<") solver for the adjoint inviscid terms."<< endl;
        if (MUSCL_AdjFlow) {
          PrintLimiterInfo(Kind_SlopeLimit_AdjFlow);
        } else {
          cout << "First order integration." << endl;
        }
      }

      cout << "The reference sharp edge distance is: " << AdjSharp_LimiterCoeff*RefElemLength*Venkat_LimiterCoeff <<". "<< endl;

    }

    if ((Kind_Solver == MAIN_SOLVER::ADJ_RANS) && (!Frozen_Visc_Cont)) {
      if (Kind_ConvNumScheme_AdjTurb == SPACE_UPWIND) {
        if (Kind_Upwind_Turb == UPWIND::SCALAR_UPWIND) cout << "Scalar upwind solver for the adjoint turbulence model." << endl;
        if (MUSCL_AdjTurb) {
          PrintLimiterInfo(Kind_SlopeLimit_AdjTurb);
        } else {
          cout << "First order integration." << endl;
        }
      }
    }

    if ((Kind_Solver == MAIN_SOLVER::NAVIER_STOKES) || (Kind_Solver == MAIN_SOLVER::RANS) ||
        (Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES) || (Kind_Solver == MAIN_SOLVER::INC_RANS) ||
        (Kind_Solver == MAIN_SOLVER::NEMO_NAVIER_STOKES) ||
        (Kind_Solver == MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_INC_RANS) ||
        (Kind_Solver == MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_RANS)) {
        cout << "Average of gradients with correction (viscous flow terms)." << endl;
    }

    if ((Kind_Solver == MAIN_SOLVER::ADJ_NAVIER_STOKES) || (Kind_Solver == MAIN_SOLVER::ADJ_RANS)) {
      cout << "Average of gradients with correction (viscous adjoint terms)." << endl;
    }

    if ((Kind_Solver == MAIN_SOLVER::RANS) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_RANS) || (Kind_Solver == MAIN_SOLVER::INC_RANS) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_INC_RANS) ) {
      cout << "Average of gradients with correction (viscous turbulence terms)." << endl;
    }

    if ((Kind_Solver == MAIN_SOLVER::ADJ_RANS) && (!Frozen_Visc_Cont)) {
      cout << "Average of gradients with correction (2nd order) for computation of adjoint viscous turbulence terms." << endl;
      if (Kind_TimeIntScheme_AdjTurb == EULER_IMPLICIT) cout << "Euler implicit method for the turbulent adjoint equation." << endl;
    }

    if(Kind_Solver != MAIN_SOLVER::FEM_EULER && Kind_Solver != MAIN_SOLVER::FEM_NAVIER_STOKES &&
       Kind_Solver != MAIN_SOLVER::FEM_RANS  && Kind_Solver != MAIN_SOLVER::FEM_LES &&
       Kind_Solver != MAIN_SOLVER::DISC_ADJ_FEM_EULER && Kind_Solver != MAIN_SOLVER::DISC_ADJ_FEM_NS &&
       Kind_Solver != MAIN_SOLVER::DISC_ADJ_FEM_RANS) {
      if (!fea){
        switch (Kind_Gradient_Method_Recon) {
          case GREEN_GAUSS: cout << "Gradient for upwind reconstruction: Green-Gauss." << endl; break;
          case LEAST_SQUARES: cout << "Gradient for upwind reconstruction: unweighted Least-Squares." << endl; break;
          case WEIGHTED_LEAST_SQUARES: cout << "Gradient for upwind reconstruction: inverse-distance weighted Least-Squares." << endl; break;
        }
        switch (Kind_Gradient_Method) {
          case GREEN_GAUSS: cout << "Gradient for viscous and source terms: Green-Gauss." << endl; break;
          case LEAST_SQUARES: cout << "Gradient for viscous and source terms: unweighted Least-Squares." << endl; break;
          case WEIGHTED_LEAST_SQUARES: cout << "Gradient for viscous and source terms: inverse-distance weighted Least-Squares." << endl; break;
        }
      }
      else{
        cout << "Spatial discretization using the Finite Element Method." << endl;
      }
    }

    if(Kind_Solver == MAIN_SOLVER::FEM_EULER || Kind_Solver == MAIN_SOLVER::FEM_NAVIER_STOKES ||
       Kind_Solver == MAIN_SOLVER::FEM_RANS  || Kind_Solver == MAIN_SOLVER::FEM_LES ||
       Kind_Solver == MAIN_SOLVER::DISC_ADJ_FEM_EULER || Kind_Solver == MAIN_SOLVER::DISC_ADJ_FEM_NS ||
       Kind_Solver == MAIN_SOLVER::DISC_ADJ_FEM_RANS) {
      if(Kind_FEM_Flow == DG) {
        cout << "Discontinuous Galerkin Finite element solver" << endl;

        switch( Riemann_Solver_FEM ) {
          case UPWIND::ROE:           cout << "Roe (with entropy fix) solver for inviscid fluxes over the faces" << endl; break;
          case UPWIND::LAX_FRIEDRICH: cout << "Lax-Friedrich solver for inviscid fluxes over the faces" << endl; break;
          case UPWIND::AUSM:          cout << "AUSM solver inviscid fluxes over the faces" << endl; break;
          case UPWIND::HLLC:          cout << "HLLC solver inviscid fluxes over the faces" << endl; break;
          default: break;
        }

        if(Kind_Solver != MAIN_SOLVER::FEM_EULER && Kind_Solver != MAIN_SOLVER::DISC_ADJ_FEM_EULER) {
          cout << "Theta symmetrizing terms interior penalty: " << Theta_Interior_Penalty_DGFEM << endl;
        }
      }

      cout << "Quadrature factor for elements with constant Jacobian:     " << Quadrature_Factor_Straight << endl;
      cout << "Quadrature factor for elements with non-constant Jacobian: " << Quadrature_Factor_Curved << endl;

      cout << "Byte alignment matrix multiplications:      " << byteAlignmentMatMul << endl;
      cout << "Padded matrix size for optimal performance: " << sizeMatMulPadding << endl;
    }

    cout << endl <<"--------------- Time Numerical Integration  ( Zone "  << iZone << " ) ------------------" << endl;

    if (!fea) {
    switch (TimeMarching) {
      case TIME_MARCHING::STEADY:
        cout << "Local time stepping (steady state simulation)." << endl; break;

      case TIME_MARCHING::TIME_STEPPING:
        cout << "Unsteady simulation using a time stepping strategy."<< endl;
        if (Unst_CFL != 0.0) {
          cout << "Time step computed by the code. Unsteady CFL number: " << Unst_CFL <<"."<< endl;
          if (Delta_UnstTime != 0.0) {
            cout << "Synchronization time provided by the user (s): "<< Delta_UnstTime << "." << endl;
          }
        }
        else cout << "Unsteady time step provided by the user (s): "<< Delta_UnstTime << "." << endl;
        break;

      case TIME_MARCHING::DT_STEPPING_1ST: case TIME_MARCHING::DT_STEPPING_2ND:
        if (TimeMarching == TIME_MARCHING::DT_STEPPING_1ST) cout << "Unsteady simulation, dual time stepping strategy (first order in time)."<< endl;
        if (TimeMarching == TIME_MARCHING::DT_STEPPING_2ND) cout << "Unsteady simulation, dual time stepping strategy (second order in time)."<< endl;
        if (Unst_CFL != 0.0) cout << "Time step computed by the code. Unsteady CFL number: " << Unst_CFL <<"."<< endl;
        else cout << "Unsteady time step provided by the user (s): "<< Delta_UnstTime << "." << endl;
        break;

      default:
        break;
      }
    }
    else {
      if (Time_Domain) {
        cout << "Dynamic structural analysis."<< endl;
        cout << "Time step provided by the user for the dynamic analysis(s): "<< Time_Step << "." << endl;
      } else {
        cout << "Static structural analysis." << endl;
      }
    }

    if ((Kind_Solver == MAIN_SOLVER::EULER) || (Kind_Solver == MAIN_SOLVER::NAVIER_STOKES) || (Kind_Solver == MAIN_SOLVER::RANS) ||
        (Kind_Solver == MAIN_SOLVER::INC_EULER) || (Kind_Solver == MAIN_SOLVER::INC_NAVIER_STOKES) || (Kind_Solver == MAIN_SOLVER::INC_RANS) ||
        (Kind_Solver == MAIN_SOLVER::NEMO_EULER) || (Kind_Solver == MAIN_SOLVER::NEMO_NAVIER_STOKES) ||
        (Kind_Solver == MAIN_SOLVER::DISC_ADJ_INC_EULER) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_INC_RANS) ||
        (Kind_Solver == MAIN_SOLVER::DISC_ADJ_EULER) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_RANS) ||
        (Kind_Solver == MAIN_SOLVER::DISC_ADJ_FEM_EULER) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_FEM_NS) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_FEM_RANS)) {
      switch (Kind_TimeIntScheme_Flow) {
        case RUNGE_KUTTA_EXPLICIT:
          cout << "Runge-Kutta explicit method for the flow equations." << endl;
          cout << "Number of steps: " << nRKStep << endl;
          cout << "Alpha coefficients: ";
          for (unsigned short iRKStep = 0; iRKStep < nRKStep; iRKStep++) {
            cout << "\t" << RK_Alpha_Step[iRKStep];
          }
          cout << endl;
          break;
        case EULER_EXPLICIT:
          cout << "Euler explicit method for the flow equations." << endl;
          break;
        case EULER_IMPLICIT:
          cout << "Euler implicit method for the flow equations." << endl;
          if (Kind_FluidModel == MUTATIONPP)
            SU2_MPI::Error("Implicit time scheme is not yet implemented with Mutation++. Use EULER_EXPLICIT.", CURRENT_FUNCTION);
          switch (Kind_Linear_Solver) {
            case BCGSTAB:
            case FGMRES:
            case RESTARTED_FGMRES:
              if (Kind_Linear_Solver == BCGSTAB)
                cout << "BCGSTAB is used for solving the linear system." << endl;
              else
                cout << "FGMRES is used for solving the linear system." << endl;
              switch (Kind_Linear_Solver_Prec) {
                case ILU: cout << "Using a ILU("<< Linear_Solver_ILU_n <<") preconditioning."<< endl; break;
                case LINELET: cout << "Using a linelet preconditioning."<< endl; break;
                case LU_SGS:  cout << "Using a LU-SGS preconditioning."<< endl; break;
                case JACOBI:  cout << "Using a Jacobi preconditioning."<< endl; break;
              }
              break;
            case SMOOTHER:
              switch (Kind_Linear_Solver_Prec) {
                case ILU:     cout << "A ILU(" << Linear_Solver_ILU_n << ")"; break;
                case LINELET: cout << "A Linelet"; break;
                case LU_SGS:  cout << "A LU-SGS"; break;
                case JACOBI:  cout << "A Jacobi"; break;
              }
              cout << " method is used for smoothing the linear system." << endl;
              break;
          }
          cout << "Convergence criteria of the linear solver: "<< Linear_Solver_Error <<"."<< endl;
          cout << "Max number of linear iterations: "<< Linear_Solver_Iter <<"."<< endl;
          break;
        case CLASSICAL_RK4_EXPLICIT:
          cout << "Classical RK4 explicit method for the flow equations." << endl;
          cout << "Number of steps: " << 4 << endl;
          cout << "Time coefficients: {0.5, 0.5, 1, 1}" << endl;
          cout << "Function coefficients: {1/6, 1/3, 1/3, 1/6}" << endl;
          break;
      }
    }

    if (fea) {
      switch (Kind_TimeIntScheme_FEA) {
        case STRUCT_TIME_INT::GENERALIZED_ALPHA:
          cout << "Generalized-alpha method." << endl;
          break;
        case STRUCT_TIME_INT::NEWMARK_IMPLICIT:
          if (Time_Domain) cout << "Newmark implicit method for the structural time integration." << endl;
          switch (Kind_Linear_Solver) {
            case BCGSTAB:
              cout << "BCGSTAB is used for solving the linear system." << endl;
              cout << "Convergence criteria of the linear solver: "<< Linear_Solver_Error <<"."<< endl;
              cout << "Max number of iterations: "<< Linear_Solver_Iter <<"."<< endl;
              break;
            case FGMRES: case RESTARTED_FGMRES:
              cout << "FGMRES is used for solving the linear system." << endl;
              cout << "Convergence criteria of the linear solver: "<< Linear_Solver_Error <<"."<< endl;
              cout << "Max number of iterations: "<< Linear_Solver_Iter <<"."<< endl;
              break;
            case CONJUGATE_GRADIENT:
              cout << "A Conjugate Gradient method is used for solving the linear system." << endl;
              cout << "Convergence criteria of the linear solver: "<< Linear_Solver_Error <<"."<< endl;
              cout << "Max number of iterations: "<< Linear_Solver_Iter <<"."<< endl;
              break;
          }
          break;
      }
    }

    if ((Kind_Solver == MAIN_SOLVER::ADJ_EULER) || (Kind_Solver == MAIN_SOLVER::ADJ_NAVIER_STOKES) || (Kind_Solver == MAIN_SOLVER::ADJ_RANS)) {
      switch (Kind_TimeIntScheme_AdjFlow) {
        case RUNGE_KUTTA_EXPLICIT:
          cout << "Runge-Kutta explicit method for the adjoint equations." << endl;
          cout << "Number of steps: " << nRKStep << endl;
          cout << "Alpha coefficients: ";
          for (unsigned short iRKStep = 0; iRKStep < nRKStep; iRKStep++) {
            cout << "\t" << RK_Alpha_Step[iRKStep];
          }
          cout << endl;
          break;
        case EULER_EXPLICIT: cout << "Euler explicit method for the adjoint equations." << endl; break;
        case EULER_IMPLICIT: cout << "Euler implicit method for the adjoint equations." << endl; break;
      }
    }

    if(Kind_Solver == MAIN_SOLVER::FEM_EULER || Kind_Solver == MAIN_SOLVER::FEM_NAVIER_STOKES ||
       Kind_Solver == MAIN_SOLVER::FEM_RANS  || Kind_Solver == MAIN_SOLVER::FEM_LES) {
      switch (Kind_TimeIntScheme_FEM_Flow) {
        case RUNGE_KUTTA_EXPLICIT:
          cout << "Runge-Kutta explicit method for the flow equations." << endl;
          cout << "Number of steps: " << nRKStep << endl;
          cout << "Alpha coefficients: ";
          for (unsigned short iRKStep = 0; iRKStep < nRKStep; iRKStep++) {
            cout << "\t" << RK_Alpha_Step[iRKStep];
          }
          cout << endl;
          break;
        case CLASSICAL_RK4_EXPLICIT:
          cout << "Classical RK4 explicit method for the flow equations." << endl;
          cout << "Number of steps: " << 4 << endl;
          cout << "Time coefficients: {0.5, 0.5, 1, 1}" << endl;
          cout << "Function coefficients: {1/6, 1/3, 1/3, 1/6}" << endl;
          break;

        case ADER_DG:
          if(nLevels_TimeAccurateLTS == 1)
            cout << "ADER-DG for the flow equations with global time stepping." << endl;
          else
            cout << "ADER-DG for the flow equations with " << nLevels_TimeAccurateLTS
                 << " levels for time accurate local time stepping." << endl;

          switch( Kind_ADER_Predictor ) {
            case ADER_ALIASED_PREDICTOR:
              cout << "An aliased approach is used in the predictor step. " << endl;
              break;
            case ADER_NON_ALIASED_PREDICTOR:
              cout << "A non-aliased approach is used in the predictor step. " << endl;
              break;
          }
          cout << "Number of time DOFs ADER-DG predictor step: " << nTimeDOFsADER_DG << endl;
          cout << "Location of time DOFs ADER-DG on the interval [-1,1]: ";
          for (unsigned short iDOF=0; iDOF<nTimeDOFsADER_DG; iDOF++) {
            cout << "\t" << TimeDOFsADER_DG[iDOF];
          }
          cout << endl;
          cout << "Time quadrature factor for ADER-DG: " << Quadrature_Factor_Time_ADER_DG << endl;
          cout << "Number of time integration points ADER-DG: " << nTimeIntegrationADER_DG << endl;
          cout << "Location of time integration points ADER-DG on the interval [-1,1]: ";
          for (unsigned short iDOF=0; iDOF<nTimeIntegrationADER_DG; iDOF++) {
            cout << "\t" << TimeIntegrationADER_DG[iDOF];
          }
          cout << endl;
          cout << "Weights of time integration points ADER-DG on the interval [-1,1]: ";
          for (unsigned short iDOF=0; iDOF<nTimeIntegrationADER_DG; iDOF++) {
            cout << "\t" << WeightsIntegrationADER_DG[iDOF];
          }
          cout << endl;
          break;
      }
    }

    if (nMGLevels !=0) {

      if (MGCycle == V_CYCLE) cout << "V Multigrid Cycle, with " << nMGLevels << " multigrid levels."<< endl;
      if (MGCycle == W_CYCLE) cout << "W Multigrid Cycle, with " << nMGLevels << " multigrid levels."<< endl;
      if (MGCycle == FULLMG_CYCLE) cout << "Full Multigrid Cycle, with " << nMGLevels << " multigrid levels."<< endl;

      cout << "Damping factor for the residual restriction: " << Damp_Res_Restric <<"."<< endl;
      cout << "Damping factor for the correction prolongation: " << Damp_Correc_Prolong <<"."<< endl;
    }

    if ((Kind_Solver != MAIN_SOLVER::FEM_ELASTICITY) && (Kind_Solver != MAIN_SOLVER::DISC_ADJ_FEM)) {

      if (!CFL_Adapt) cout << "No CFL adaptation." << endl;
      else cout << "CFL adaptation. Factor down: "<< CFL_AdaptParam[0] <<", factor up: "<< CFL_AdaptParam[1]
        <<",\n                lower limit: "<< CFL_AdaptParam[2] <<", upper limit: " << CFL_AdaptParam[3]
        <<",\n                acceptable linear residual: "<< CFL_AdaptParam[4]
        <<"'\n                starting iteration: "<< CFL_AdaptParam[5] << "." << endl;

      if (nMGLevels !=0) {
        PrintingToolbox::CTablePrinter MGTable(&std::cout);

        MGTable.AddColumn("MG Level",         10);
        MGTable.AddColumn("Presmooth",     10);
        MGTable.AddColumn("PostSmooth",    10);
        MGTable.AddColumn("CorrectSmooth", 10);
        MGTable.SetAlign(PrintingToolbox::CTablePrinter::RIGHT);
        MGTable.PrintHeader();
        for (unsigned short iLevel = 0; iLevel < nMGLevels+1; iLevel++) {
          MGTable << iLevel << MG_PreSmooth[iLevel] << MG_PostSmooth[iLevel] << MG_CorrecSmooth[iLevel];
        }
        MGTable.PrintFooter();
      }
      if (TimeMarching != TIME_MARCHING::TIME_STEPPING) {
        cout << "Courant-Friedrichs-Lewy number:   ";
        cout.precision(3);
        cout.width(6); cout << CFL[0];
        cout << endl;
      }

    }

    if ((Kind_Solver == MAIN_SOLVER::RANS) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_RANS) ||
        (Kind_Solver == MAIN_SOLVER::INC_RANS) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_INC_RANS))
      if (Kind_TimeIntScheme_Turb == EULER_IMPLICIT)
        cout << "Euler implicit time integration for the turbulence model." << endl;
  }

  if (val_software == SU2_COMPONENT::SU2_CFD) {

    cout << endl <<"------------------ Convergence Criteria  ( Zone "  << iZone << " ) ---------------------" << endl;

    cout << "Maximum number of solver subiterations: " << nInnerIter <<"."<< endl;
    if (Multizone_Problem)
      cout << "Maximum number of solver outer iterations: " << nOuterIter <<"."<< endl;
    if (Time_Domain)
      cout << "Maximum number of physical time-steps: " << nTimeIter <<"."<< endl;

    cout << "Begin convergence monitoring at iteration " << StartConv_Iter << "." << endl;
    cout << "Residual minimum value: 1e" << MinLogResidual << "." << endl;
    cout << "Cauchy series min. value: " << Cauchy_Eps << "." << endl;
    cout << "Number of Cauchy elements: " << Cauchy_Elems << "." << endl;
    if(Cauchy_Elems <1){
      SU2_MPI::Error(to_string(Cauchy_Elems) + string(" Cauchy elements are no viable input. Please check your configuration file."), CURRENT_FUNCTION);
    }
    cout << "Begin windowed time average at iteration " << StartWindowIteration << "." << endl;

    if(Wnd_Cauchy_Crit){
      cout << "Begin time convergence monitoring at iteration " << Wnd_StartConv_Iter + StartWindowIteration << "." << endl;
      cout << "Time cauchy series min. value: " << Wnd_Cauchy_Eps << "." << endl;
      cout << "Number of Cauchy elements: " << Wnd_Cauchy_Elems << "." << endl;
      if(Wnd_Cauchy_Elems <1){
        SU2_MPI::Error(to_string(Wnd_Cauchy_Elems) +string(" Cauchy elements are no viable input. Please check your configuration file."), CURRENT_FUNCTION);
      }
    }
  }

  cout << endl <<"-------------------- Output Information ( Zone "  << iZone << " ) ----------------------" << endl;

  if (val_software == SU2_COMPONENT::SU2_CFD) {

    if (nVolumeOutputFiles != 0) {
      cout << "File writing frequency: " << endl;
      PrintingToolbox::CTablePrinter FileFreqTable(&std::cout);
      FileFreqTable.AddColumn("File", 25);
      FileFreqTable.AddColumn("Frequency", 10);
      FileFreqTable.SetAlign(PrintingToolbox::CTablePrinter::RIGHT);
      FileFreqTable.PrintHeader();

      for (auto iFreq = 0; iFreq < nVolumeOutputFiles; iFreq++){
        /*--- find the key belonging to the value in the map---*/
        for (auto& it : Output_Map) {
          if (it.second == VolumeOutputFiles[iFreq]) {
            FileFreqTable << it.first << VolumeOutputFrequencies[iFreq];
            break;
          }
        }
      }

      FileFreqTable.PrintFooter();
    }

    cout << "Writing the convergence history file every " << HistoryWrtFreq[2] <<" inner iterations."<< endl;
    if (Multizone_Problem){
      cout << "Writing the convergence history file every " << HistoryWrtFreq[1] <<" outer iterations."<< endl;
    }
    if (Time_Domain) {
      cout << "Writing the convergence history file every " << HistoryWrtFreq[0] <<" time iterations."<< endl;
    }
    cout << "Writing the screen convergence history every " << ScreenWrtFreq[2] <<" inner iterations."<< endl;
    if (Multizone_Problem){
      cout << "Writing the screen convergence history every " << ScreenWrtFreq[1] <<" outer iterations."<< endl;
    }
    if (Time_Domain) {
      cout << "Writing the screen convergence history every " << ScreenWrtFreq[0] <<" time iterations."<< endl;
    }

    switch (Tab_FileFormat) {
      case TAB_OUTPUT::TAB_CSV: cout << "The tabular file format is CSV (.csv)." << endl; break;
      case TAB_OUTPUT::TAB_TECPLOT: cout << "The tabular file format is Tecplot (.dat)." << endl; break;
    }

    cout << "Convergence history file name: " << Conv_FileName << "." << endl;

    cout << "Forces breakdown file name: " << Breakdown_FileName << "." << endl;


    if (!ContinuousAdjoint && !DiscreteAdjoint) {
      cout << "Surface file name: " << SurfCoeff_FileName << "." << endl;
      cout << "Volume file name: " << Volume_FileName << "." << endl;
      cout << "Restart file name: " << Restart_FileName << "." << endl;
    }

    if (ContinuousAdjoint || DiscreteAdjoint) {
      cout << "Adjoint solution file name: " << Solution_AdjFileName << "." << endl;
      cout << "Restart adjoint file name: " << Restart_AdjFileName << "." << endl;
      cout << "Adjoint variables file name: " << Adj_FileName << "." << endl;
      cout << "Surface adjoint file name: " << SurfAdjCoeff_FileName << "." << endl;
    }

  }

  if (val_software == SU2_COMPONENT::SU2_SOL) {
    switch (Tab_FileFormat) {
      case TAB_OUTPUT::TAB_CSV: cout << "The tabular file format is CSV (.csv)." << endl; break;
      case TAB_OUTPUT::TAB_TECPLOT: cout << "The tabular file format is Tecplot (.dat)." << endl; break;
    }
    cout << "Flow variables file name: " << Volume_FileName << "." << endl;
  }

  if (val_software == SU2_COMPONENT::SU2_DEF) {
    cout << "Output mesh file name: " << Mesh_Out_FileName << ". " << endl;
    switch (GetDeform_Stiffness_Type()) {
      case INVERSE_VOLUME:
        cout << "Cell stiffness scaled by inverse of the cell volume." << endl;
        break;
      case SOLID_WALL_DISTANCE:
        cout << "Cell stiffness scaled by distance to nearest solid surface." << endl;
        break;
      case CONSTANT_STIFFNESS:
        cout << "Imposing constant cell stiffness." << endl;
        break;
    }
  }

  if (val_software == SU2_COMPONENT::SU2_DOT) {
    if (DiscreteAdjoint) {
      cout << "Output Volume Sensitivity file name: " << VolSens_FileName << ". " << endl;
      cout << "Output Surface Sensitivity file name: " << SurfSens_FileName << ". " << endl;
    }
    cout << "Output gradient file name: " << ObjFunc_Grad_FileName << ". " << endl;
  }

  cout << endl <<"------------- Config File Boundary Information ( Zone "  << iZone << " ) ---------------" << endl;

  PrintingToolbox::CTablePrinter BoundaryTable(&std::cout);
  BoundaryTable.AddColumn("Marker Type", 35);
  BoundaryTable.AddColumn("Marker Name", 35);

  BoundaryTable.PrintHeader();

  if (nMarker_Euler != 0) {
    BoundaryTable << "Euler wall";
    for (iMarker_Euler = 0; iMarker_Euler < nMarker_Euler; iMarker_Euler++) {
      BoundaryTable << Marker_Euler[iMarker_Euler];
      if (iMarker_Euler < nMarker_Euler-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_FarField != 0) {
    BoundaryTable << "Far-field";
    for (iMarker_FarField = 0; iMarker_FarField < nMarker_FarField; iMarker_FarField++) {
      BoundaryTable << Marker_FarField[iMarker_FarField];
      if (iMarker_FarField < nMarker_FarField-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_SymWall != 0) {
    BoundaryTable << "Symmetry plane";
    for (iMarker_SymWall = 0; iMarker_SymWall < nMarker_SymWall; iMarker_SymWall++) {
      BoundaryTable << Marker_SymWall[iMarker_SymWall];
      if (iMarker_SymWall < nMarker_SymWall-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_PerBound != 0) {
    BoundaryTable << "Periodic boundary";
    for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++) {
      BoundaryTable << Marker_PerBound[iMarker_PerBound];
      if (iMarker_PerBound < nMarker_PerBound-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_NearFieldBound != 0) {
    BoundaryTable << "Near-field boundary";
    for (iMarker_NearFieldBound = 0; iMarker_NearFieldBound < nMarker_NearFieldBound; iMarker_NearFieldBound++) {
      BoundaryTable << Marker_NearFieldBound[iMarker_NearFieldBound];
      if (iMarker_NearFieldBound < nMarker_NearFieldBound-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Deform_Mesh != 0) {
    BoundaryTable << "Deformable mesh boundary";
    for (iMarker_Deform_Mesh = 0; iMarker_Deform_Mesh < nMarker_Deform_Mesh; iMarker_Deform_Mesh++) {
      BoundaryTable << Marker_Deform_Mesh[iMarker_Deform_Mesh];
      if (iMarker_Deform_Mesh < nMarker_Deform_Mesh-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Deform_Mesh_Sym_Plane != 0) {
    BoundaryTable << "Symmetric deformable mesh boundary";
    for (iMarker_Deform_Mesh_Sym_Plane = 0; iMarker_Deform_Mesh_Sym_Plane < nMarker_Deform_Mesh_Sym_Plane; iMarker_Deform_Mesh_Sym_Plane++) {
      BoundaryTable << Marker_Deform_Mesh_Sym_Plane[iMarker_Deform_Mesh_Sym_Plane];
      if (iMarker_Deform_Mesh_Sym_Plane < nMarker_Deform_Mesh_Sym_Plane-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Fluid_Load != 0) {
    BoundaryTable << "Fluid loads boundary";
    for (iMarker_Fluid_Load = 0; iMarker_Fluid_Load < nMarker_Fluid_Load; iMarker_Fluid_Load++) {
      BoundaryTable << Marker_Fluid_Load[iMarker_Fluid_Load];
      if (iMarker_Fluid_Load < nMarker_Fluid_Load-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Fluid_InterfaceBound != 0) {
    BoundaryTable << "Fluid interface boundary";
    for (iMarker_Fluid_InterfaceBound = 0; iMarker_Fluid_InterfaceBound < nMarker_Fluid_InterfaceBound; iMarker_Fluid_InterfaceBound++) {
      BoundaryTable << Marker_Fluid_InterfaceBound[iMarker_Fluid_InterfaceBound];
      if (iMarker_Fluid_InterfaceBound < nMarker_Fluid_InterfaceBound-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Internal != 0) {
    BoundaryTable << "Internal boundary";
    for (iMarker_Internal = 0; iMarker_Internal < nMarker_Internal; iMarker_Internal++) {
      BoundaryTable << Marker_Internal[iMarker_Internal];
      if (iMarker_Internal < nMarker_Internal-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Inlet != 0) {
    BoundaryTable << "Inlet boundary";
    for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++) {
      BoundaryTable << Marker_Inlet[iMarker_Inlet];
      if (iMarker_Inlet < nMarker_Inlet-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Inlet_Species != 0) {
    BoundaryTable << "Species Inlet boundary";
    for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet_Species; iMarker_Inlet++) {
      BoundaryTable << Marker_Inlet_Species[iMarker_Inlet];
      if (iMarker_Inlet < nMarker_Inlet_Species-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Riemann != 0) {
    BoundaryTable << "Riemann boundary";
    for (iMarker_Riemann = 0; iMarker_Riemann < nMarker_Riemann; iMarker_Riemann++) {
      BoundaryTable << Marker_Riemann[iMarker_Riemann];
      if (iMarker_Riemann < nMarker_Riemann-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Giles != 0) {
    BoundaryTable << "Giles boundary";
    for (iMarker_Giles = 0; iMarker_Giles < nMarker_Giles; iMarker_Giles++) {
      BoundaryTable << Marker_Giles[iMarker_Giles];
      if (iMarker_Giles < nMarker_Giles-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_MixingPlaneInterface != 0) {
    BoundaryTable << "MixingPlane boundary";
    for (iMarker_MixingPlaneInterface = 0; iMarker_MixingPlaneInterface < nMarker_MixingPlaneInterface; iMarker_MixingPlaneInterface++) {
      BoundaryTable << Marker_MixingPlaneInterface[iMarker_MixingPlaneInterface];
      if (iMarker_MixingPlaneInterface < nMarker_MixingPlaneInterface-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_EngineInflow != 0) {
    BoundaryTable << "Engine inflow boundary";
    for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++) {
      BoundaryTable << Marker_EngineInflow[iMarker_EngineInflow];
      if (iMarker_EngineInflow < nMarker_EngineInflow-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_EngineExhaust != 0) {
    BoundaryTable << "Engine exhaust boundary";
    for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++) {
      BoundaryTable << Marker_EngineExhaust[iMarker_EngineExhaust];
      if (iMarker_EngineExhaust < nMarker_EngineExhaust-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Supersonic_Inlet != 0) {
    BoundaryTable << "Supersonic inlet boundary";
    for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++) {
      BoundaryTable << Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet];
      if (iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Supersonic_Outlet != 0) {
    BoundaryTable << "Supersonic outlet boundary";
    for (iMarker_Supersonic_Outlet = 0; iMarker_Supersonic_Outlet < nMarker_Supersonic_Outlet; iMarker_Supersonic_Outlet++) {
      BoundaryTable << Marker_Supersonic_Outlet[iMarker_Supersonic_Outlet];
      if (iMarker_Supersonic_Outlet < nMarker_Supersonic_Outlet-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Outlet != 0) {
    BoundaryTable << "Outlet boundary";
    for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++) {
      BoundaryTable << Marker_Outlet[iMarker_Outlet];
      if (iMarker_Outlet < nMarker_Outlet-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Isothermal != 0) {
    BoundaryTable << "Isothermal wall";
    for (iMarker_Isothermal = 0; iMarker_Isothermal < nMarker_Isothermal; iMarker_Isothermal++) {
      BoundaryTable << Marker_Isothermal[iMarker_Isothermal];
      if (iMarker_Isothermal < nMarker_Isothermal-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Smoluchowski_Maxwell != 0) {
    BoundaryTable << "Smoluchowski/Maxwell jump wall";
    for (iMarker_Smoluchowski_Maxwell = 0; iMarker_Smoluchowski_Maxwell < nMarker_Smoluchowski_Maxwell; iMarker_Smoluchowski_Maxwell++) {
      BoundaryTable << Marker_Smoluchowski_Maxwell[iMarker_Smoluchowski_Maxwell];
      if (iMarker_Smoluchowski_Maxwell< nMarker_Smoluchowski_Maxwell-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_HeatFlux != 0) {
    BoundaryTable << "Heat flux wall";
    for (iMarker_HeatFlux = 0; iMarker_HeatFlux < nMarker_HeatFlux; iMarker_HeatFlux++) {
      BoundaryTable << Marker_HeatFlux[iMarker_HeatFlux];
      if (iMarker_HeatFlux < nMarker_HeatFlux-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_HeatTransfer != 0) {
    BoundaryTable << "Heat transfer wall";
    for (iMarker_HeatTransfer = 0; iMarker_HeatTransfer < nMarker_HeatTransfer; iMarker_HeatTransfer++) {
      BoundaryTable << Marker_HeatTransfer[iMarker_HeatTransfer];
      if (iMarker_HeatTransfer < nMarker_HeatTransfer-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nWall_Catalytic != 0) {
    BoundaryTable << "Catalytic wall";
    for (iWall_Catalytic = 0; iWall_Catalytic < nWall_Catalytic; iWall_Catalytic++) {
      BoundaryTable << Wall_Catalytic[iWall_Catalytic];
      if (iWall_Catalytic < nWall_Catalytic-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Clamped != 0) {
    BoundaryTable << "Clamped boundary";
    for (iMarker_Clamped = 0; iMarker_Clamped < nMarker_Clamped; iMarker_Clamped++) {
      BoundaryTable << Marker_Clamped[iMarker_Clamped];
      if (iMarker_Clamped < nMarker_Clamped-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Displacement != 0) {
    BoundaryTable << "Displacement boundary";
    for (iMarker_Displacement = 0; iMarker_Displacement < nMarker_Displacement; iMarker_Displacement++) {
      BoundaryTable << Marker_Displacement[iMarker_Displacement];
      if (iMarker_Displacement < nMarker_Displacement-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Load != 0) {
    BoundaryTable << "Normal load boundary";
    for (iMarker_Load = 0; iMarker_Load < nMarker_Load; iMarker_Load++) {
      BoundaryTable << Marker_Load[iMarker_Load];
      if (iMarker_Load < nMarker_Load-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Damper != 0) {
    BoundaryTable << "Damper boundary";
    for (iMarker_Damper = 0; iMarker_Damper < nMarker_Damper; iMarker_Damper++) {
      BoundaryTable << Marker_Damper[iMarker_Damper];
      if (iMarker_Damper < nMarker_Damper-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Load_Dir != 0) {
    BoundaryTable << "Load boundary";
    for (iMarker_Load_Dir = 0; iMarker_Load_Dir < nMarker_Load_Dir; iMarker_Load_Dir++) {
      BoundaryTable << Marker_Load_Dir[iMarker_Load_Dir];
      if (iMarker_Load_Dir < nMarker_Load_Dir-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Disp_Dir != 0) {
    BoundaryTable << "Disp boundary";
    for (iMarker_Disp_Dir = 0; iMarker_Disp_Dir < nMarker_Disp_Dir; iMarker_Disp_Dir++) {
      BoundaryTable << Marker_Disp_Dir[iMarker_Disp_Dir];
      if (iMarker_Disp_Dir < nMarker_Disp_Dir-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Emissivity != 0) {
    BoundaryTable << "Radiative boundary";
    for (iMarker_Emissivity = 0; iMarker_Emissivity < nMarker_Emissivity; iMarker_Emissivity++) {
      BoundaryTable << Marker_Emissivity[iMarker_Emissivity]; // << "(" << Wall_Emissivity[iMarker_Emissivity] << ")";
      if (iMarker_Emissivity < nMarker_Emissivity-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_StrongBC != 0) {
    BoundaryTable << "Strong boundary";
    for (iMarker_StrongBC = 0; iMarker_StrongBC < nMarker_StrongBC; iMarker_StrongBC++) {
      BoundaryTable << Marker_StrongBC[iMarker_StrongBC];
      if (iMarker_StrongBC < nMarker_StrongBC-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_Custom != 0) {
    BoundaryTable << "Custom boundary";
    for (iMarker_Custom = 0; iMarker_Custom < nMarker_Custom; iMarker_Custom++) {
      BoundaryTable << Marker_Custom[iMarker_Custom];
      if (iMarker_Custom < nMarker_Custom-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_ActDiskInlet != 0) {
    BoundaryTable << "Actuator disk (inlet) boundary";
    for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++) {
      BoundaryTable << Marker_ActDiskInlet[iMarker_ActDiskInlet];
      if (iMarker_ActDiskInlet < nMarker_ActDiskInlet-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_ActDiskOutlet != 0) {
    BoundaryTable << "Actuator disk (outlet) boundary";
    for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++) {
      BoundaryTable << Marker_ActDiskOutlet[iMarker_ActDiskOutlet];
      if (iMarker_ActDiskOutlet < nMarker_ActDiskOutlet-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_SobolevBC != 0) {
    BoundaryTable << "Sobolev boundary";
    for (iMarker_SobolevBC = 0; iMarker_SobolevBC < nMarker_SobolevBC; iMarker_SobolevBC++) {
      BoundaryTable << Marker_SobolevBC[iMarker_SobolevBC];
      if (iMarker_SobolevBC < nMarker_SobolevBC-1)  BoundaryTable << " ";
    }
    BoundaryTable.PrintFooter();
  }

  if (nMarker_ActDiskOutlet != 0) {
    if (GetKind_ActDisk() == VARIABLE_LOAD) {
      cout << endl << "Actuator disk with variable load." << endl;
      cout << "Actuator disk data read from file: " << GetActDisk_FileName() << endl;
    }
  }

  if (nMarker_ActDiskOutlet != 0) {
    if (GetKind_ActDisk() == BLADE_ELEMENT) {
      cout << endl << "Actuator disk with blade element momentum (BEM) method." << endl;
      cout << "Actuator disk BEM method propeller data read from file: " << GetBEM_prop_filename() << endl;
    }
  }
}

bool CConfig::TokenizeString(string & str, string & option_name,
                             vector<string> & option_value) {
  const string delimiters(" (){}:,\t\n\v\f\r");
  // check for comments or empty string
  string::size_type pos, last_pos;
  pos = str.find_first_of('%');
  if ( (str.length() == 0) || (pos == 0) ) {
    // str is empty or a comment line, so no option here
    return false;
  }
  if (pos != string::npos) {
    // remove comment at end if necessary
    str.erase(pos);
  }

  // look for line composed on only delimiters (usually whitespace)
  pos = str.find_first_not_of(delimiters);
  if (pos == string::npos) {
    return false;
  }

  // find the equals sign and split string
  string name_part, value_part;
  pos = str.find('=');
  if (pos == string::npos) {
    cerr << "Error in TokenizeString(): "
    << "line in the configuration file with no \"=\" sign."
    << endl;
    cout << "Look for: " << str << endl;
    cout << "str.length() = " << str.length() << endl;
    throw(-1);
  }
  name_part = str.substr(0, pos);
  value_part = str.substr(pos+1, string::npos);
  //cout << "name_part  = |" << name_part  << "|" << endl;
  //cout << "value_part = |" << value_part << "|" << endl;

  // the first_part should consist of one string with no interior delimiters
  last_pos = name_part.find_first_not_of(delimiters, 0);
  pos = name_part.find_first_of(delimiters, last_pos);
  if ( (name_part.length() == 0) || (last_pos == string::npos) ) {
    cerr << "Error in CConfig::TokenizeString(): "
    << "line in the configuration file with no name before the \"=\" sign."
    << endl;
    throw(-1);
  }
  if (pos == string::npos) pos = name_part.length();
  option_name = name_part.substr(last_pos, pos - last_pos);
  last_pos = name_part.find_first_not_of(delimiters, pos);
  if (last_pos != string::npos) {
    cerr << "Error in TokenizeString(): "
    << "two or more options before an \"=\" sign in the configuration file."
    << endl;
    throw(-1);
  }
  StringToUpperCase(option_name);

  //cout << "option_name = |" << option_name << "|" << endl;
  //cout << "pos = " << pos << ": last_pos = " << last_pos << endl;

  // now fill the option value vector
  option_value.clear();
  last_pos = value_part.find_first_not_of(delimiters, 0);

  // detect a raw string
  if (value_part[last_pos] == '\'' && value_part.back() == '\'') {
    option_value.push_back(value_part.substr(last_pos+1, value_part.size()-last_pos-2));
    return true;
  }

  pos = value_part.find_first_of(delimiters, last_pos);
  while (string::npos != pos || string::npos != last_pos) {
    // add token to the vector<string>
    option_value.push_back(value_part.substr(last_pos, pos - last_pos));
    // skip delimiters
    last_pos = value_part.find_first_not_of(delimiters, pos);
    // find next "non-delimiter"
    pos = value_part.find_first_of(delimiters, last_pos);
  }
  if (option_value.empty()) {
    cerr << "Error in TokenizeString(): "
    << "option " << option_name << " in configuration file with no value assigned."
    << endl;
    throw(-1);
  }

#if 0
  cout << "option value(s) = ";
  for (unsigned int i = 0; i < option_value.size(); i++)
    cout << option_value[i] << " ";
  cout << endl;
#endif

  // look for ';' DV delimiters attached to values
  vector<string>::iterator it;
  it = option_value.begin();
  while (it != option_value.end()) {
    if (it->compare(";") == 0) {
      it++;
      continue;
    }

    pos = it->find(';');
    if (pos != string::npos) {
      string before_semi = it->substr(0, pos);
      string after_semi= it->substr(pos+1, string::npos);
      if (before_semi.empty()) {
        *it = ";";
        it++;
        option_value.insert(it, after_semi);
      } else {
        *it = before_semi;
        it++;
        vector<string> to_insert;
        to_insert.emplace_back(";");
        if (!after_semi.empty())
          to_insert.push_back(after_semi);
        option_value.insert(it, to_insert.begin(), to_insert.end());
      }
      it = option_value.begin(); // go back to beginning; not efficient
      continue;
    }       it++;

  }
#if 0
  cout << "option value(s) = ";
  for (unsigned int i = 0; i < option_value.size(); i++)
    cout << option_value[i] << " ";
  cout << endl;
#endif
  // remove any consecutive ";"
  it = option_value.begin();
  bool semi_at_prev = false;
  while (it != option_value.end()) {
    if (semi_at_prev) {
      if (it->compare(";") == 0) {
        option_value.erase(it);
        it = option_value.begin();
        semi_at_prev = false;
        continue;
      }
    }
    semi_at_prev = it->compare(";") == 0;
    it++;
  }

#if 0
  cout << "option value(s) = ";
  for (unsigned int i = 0; i < option_value.size(); i++)
    cout << option_value[i] << " ";
  cout << endl;
#endif
  return true;
}

unsigned short CConfig::GetMarker_CfgFile_TagBound(const string& val_marker) const {

  unsigned short iMarker_CfgFile;

  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++) {
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker)
      return iMarker_CfgFile;
  }
  SU2_MPI::Error(string("The configuration file doesn't have any definition for marker ") + val_marker, CURRENT_FUNCTION);
  return 0;
}

string CConfig::GetMarker_CfgFile_TagBound(unsigned short val_marker) const {
  return Marker_CfgFile_TagBound[val_marker];
}

unsigned short CConfig::GetMarker_CfgFile_KindBC(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_KindBC[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_Monitoring(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_Monitoring[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_GeoEval(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_GeoEval[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_Designing(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_Designing[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_Plotting(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_Plotting[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_Analyze(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_Analyze[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_ZoneInterface(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_ZoneInterface[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_Turbomachinery(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_Turbomachinery[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_TurbomachineryFlag(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_TurbomachineryFlag[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_MixingPlaneInterface(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_MixingPlaneInterface[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_Giles(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_Giles[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_DV(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_DV[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_Moving(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_Moving[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_Deform_Mesh(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_Deform_Mesh[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_Deform_Mesh_Sym_Plane(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_Deform_Mesh_Sym_Plane[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_Fluid_Load(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_Fluid_Load[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_PyCustom(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile=0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_PyCustom[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_PerBound(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_PerBound[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_ZoneInterface(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_ZoneInterface[iMarker_CfgFile];
}

unsigned short CConfig::GetMarker_CfgFile_SobolevBC(const string& val_marker) const {
  unsigned short iMarker_CfgFile;
  for (iMarker_CfgFile = 0; iMarker_CfgFile < nMarker_CfgFile; iMarker_CfgFile++)
    if (Marker_CfgFile_TagBound[iMarker_CfgFile] == val_marker) break;
  return Marker_CfgFile_SobolevBC[iMarker_CfgFile];
}

bool CConfig::GetViscous_Wall(unsigned short iMarker) const {

  return (Marker_All_KindBC[iMarker] == HEAT_FLUX  ||
          Marker_All_KindBC[iMarker] == ISOTHERMAL ||
          Marker_All_KindBC[iMarker] == HEAT_TRANSFER ||
          Marker_All_KindBC[iMarker] == SMOLUCHOWSKI_MAXWELL ||
          Marker_All_KindBC[iMarker] == CHT_WALL_INTERFACE);
}

bool CConfig::GetCatalytic_Wall(unsigned short iMarker) const {

  bool catalytic = false;
  for (unsigned short iMarker_Catalytic = 0; iMarker_Catalytic < nWall_Catalytic; iMarker_Catalytic++){
    string Catalytic_Tag = Wall_Catalytic[iMarker_Catalytic];
    if (Catalytic_Tag == Marker_All_TagBound[iMarker]) { catalytic = true; break; }
  }

  return catalytic;
}

bool CConfig::GetSolid_Wall(unsigned short iMarker) const {

  return GetViscous_Wall(iMarker) ||
         Marker_All_KindBC[iMarker] == EULER_WALL;
}

void CConfig::SetSurface_Movement(unsigned short iMarker, unsigned short kind_movement) {

  auto* new_surface_movement = new unsigned short[nMarker_Moving + 1];
  auto* new_marker_moving = new string[nMarker_Moving+1];

  for (unsigned short iMarker_Moving = 0; iMarker_Moving < nMarker_Moving; iMarker_Moving++){
    new_surface_movement[iMarker_Moving] = Kind_SurfaceMovement[iMarker_Moving];
    new_marker_moving[iMarker_Moving] = Marker_Moving[iMarker_Moving];
  }

  if (nKind_SurfaceMovement > 0){
    delete [] Marker_Moving;
    delete [] Kind_SurfaceMovement;
  }

  Kind_SurfaceMovement = new_surface_movement;
  Marker_Moving        = new_marker_moving;

  Kind_SurfaceMovement[nMarker_Moving] = kind_movement;
  Marker_Moving[nMarker_Moving] = Marker_All_TagBound[iMarker];

  nMarker_Moving++;
  nKind_SurfaceMovement++;

}

CConfig::~CConfig() {

  unsigned long iDV, iMarker;

  /*--- Delete all of the option objects in the global option map ---*/

  for(auto itr = option_map.begin(); itr != option_map.end(); itr++) {
    delete itr->second;
  }

  delete [] TimeDOFsADER_DG;
  delete [] TimeIntegrationADER_DG;
  delete [] WeightsIntegrationADER_DG;

  /*--- Free memory for Aeroelastic problems. ---*/

  delete[] Aeroelastic_pitch;
  delete[] Aeroelastic_plunge;

  /*--- Marker pointers ---*/

  delete[] Marker_CfgFile_GeoEval;
  delete[] Marker_All_GeoEval;

  delete[] Marker_CfgFile_TagBound;
  delete[] Marker_All_TagBound;

  delete[] Marker_CfgFile_KindBC;
  delete[] Marker_All_KindBC;

  delete[] Marker_CfgFile_Monitoring;
  delete[] Marker_All_Monitoring;

  delete[] Marker_CfgFile_Designing;
  delete[] Marker_All_Designing;

  delete[] Marker_CfgFile_Plotting;
  delete[] Marker_All_Plotting;

  delete[] Marker_CfgFile_Analyze;
  delete[] Marker_All_Analyze;

  delete[] Marker_CfgFile_ZoneInterface;
  delete[] Marker_All_ZoneInterface;

  delete[] Marker_CfgFile_DV;
  delete[] Marker_All_DV;

  delete[] Marker_CfgFile_Moving;
  delete[] Marker_All_Moving;

  delete[] Marker_CfgFile_Deform_Mesh;
  delete[] Marker_All_Deform_Mesh;

  delete[] Marker_CfgFile_Deform_Mesh_Sym_Plane;
  delete[] Marker_All_Deform_Mesh_Sym_Plane;

  delete[] Marker_CfgFile_Fluid_Load;
  delete[] Marker_All_Fluid_Load;

  delete[] Marker_CfgFile_PyCustom;
  delete[] Marker_All_PyCustom;

  delete[] Marker_CfgFile_PerBound;
  delete[] Marker_All_PerBound;

  delete[] Marker_CfgFile_Turbomachinery;
  delete[] Marker_All_Turbomachinery;

  delete[] Marker_CfgFile_TurbomachineryFlag;
  delete[] Marker_All_TurbomachineryFlag;

  delete[] Marker_CfgFile_Giles;
  delete[] Marker_All_Giles;

  delete[] Marker_CfgFile_MixingPlaneInterface;
  delete[] Marker_All_MixingPlaneInterface;

  delete[] Marker_CfgFile_SobolevBC;
  delete[] Marker_All_SobolevBC;

  delete[] Marker_All_SendRecv;

  if (DV_Value != nullptr) {
    for (iDV = 0; iDV < nDV; iDV++) delete[] DV_Value[iDV];
    delete [] DV_Value;
  }

  if (ParamDV != nullptr) {
    for (iDV = 0; iDV < nDV; iDV++) delete[] ParamDV[iDV];
    delete [] ParamDV;
  }

  delete[] Exhaust_Pressure;
  delete[] Exhaust_Temperature;
  delete[] Exhaust_MassFlow;
  delete[] Exhaust_TotalPressure;
  delete[] Exhaust_TotalTemperature;
  delete[] Exhaust_GrossThrust;
  delete[] Exhaust_Force;
  delete[] Exhaust_Power;

  delete[] Inflow_Mach;
  delete[] Inflow_Pressure;
  delete[] Inflow_MassFlow;
  delete[] Inflow_ReverseMassFlow;
  delete[] Inflow_TotalPressure;
  delete[] Inflow_Temperature;
  delete[] Inflow_TotalTemperature;
  delete[] Inflow_RamDrag;
  delete[] Inflow_Force;
  delete[] Inflow_Power;

  delete[] Engine_Power;
  delete[] Engine_Mach;
  delete[] Engine_Force;
  delete[] Engine_NetThrust;
  delete[] Engine_GrossThrust;
  delete[] Engine_Area;

  delete[] ActDiskInlet_MassFlow;
  delete[] ActDiskInlet_Temperature;
  delete[] ActDiskInlet_TotalTemperature;
  delete[] ActDiskInlet_Pressure;
  delete[] ActDiskInlet_TotalPressure;
  delete[] ActDiskInlet_RamDrag;
  delete[] ActDiskInlet_Force;
  delete[] ActDiskInlet_Power;

  delete[] ActDiskOutlet_MassFlow;
  delete[] ActDiskOutlet_Temperature;
  delete[] ActDiskOutlet_TotalTemperature;
  delete[] ActDiskOutlet_Pressure;
  delete[] ActDiskOutlet_TotalPressure;
  delete[] ActDiskOutlet_GrossThrust;
  delete[] ActDiskOutlet_Force;
  delete[] ActDiskOutlet_Power;

  delete[] ActDiskOutlet_Thrust_BEM;
  delete[] ActDiskOutlet_Torque_BEM;

  delete[] Outlet_MassFlow;
  delete[] Outlet_Density;
  delete[] Outlet_Area;

  delete[] ActDisk_DeltaPress;
  delete[] ActDisk_DeltaTemp;
  delete[] ActDisk_TotalPressRatio;
  delete[] ActDisk_TotalTempRatio;
  delete[] ActDisk_StaticPressRatio;
  delete[] ActDisk_StaticTempRatio;
  delete[] ActDisk_Power;
  delete[] ActDisk_MassFlow;
  delete[] ActDisk_Mach;
  delete[] ActDisk_Force;
  delete[] ActDisk_NetThrust;
  delete[] ActDisk_BCThrust;
  delete[] ActDisk_BCThrust_Old;
  delete[] ActDisk_GrossThrust;
  delete[] ActDisk_Area;
  delete[] ActDisk_ReverseMassFlow;

  delete[] Surface_MassFlow;
  delete[] Surface_Mach;
  delete[] Surface_Temperature;
  delete[] Surface_Pressure;
  delete[] Surface_Density;
  delete[] Surface_Enthalpy;
  delete[] Surface_NormalVelocity;
  delete[] Surface_Uniformity;
  delete[] Surface_SecondaryStrength;
  delete[] Surface_SecondOverUniform;
  delete[] Surface_MomentumDistortion;
  delete[] Surface_TotalTemperature;
  delete[] Surface_TotalPressure;
  delete[] Surface_PressureDrop;
  delete[] Surface_Species_0;
  delete[] Surface_Species_Variance;
  delete[] Surface_DC60;
  delete[] Surface_IDC;
  delete[] Surface_IDC_Mach;
  delete[] Surface_IDR;

  if (Riemann_FlowDir != nullptr) {
    for (iMarker = 0; iMarker < nMarker_Riemann; iMarker++)
      delete [] Riemann_FlowDir[iMarker];
    delete [] Riemann_FlowDir;
  }

  if (Giles_FlowDir != nullptr) {
    for (iMarker = 0; iMarker < nMarker_Giles; iMarker++)
      delete [] Giles_FlowDir[iMarker];
    delete [] Giles_FlowDir;
  }

  delete[] PlaneTag;
  delete[] CFL;

  /*--- Delete some arrays needed just for initializing options. ---*/

  delete [] FFDTag;
  delete [] nDV_Value;

  delete [] Kind_Data_Riemann;
  delete [] Riemann_Var1;
  delete [] Riemann_Var2;
  delete [] Kind_Data_Giles;
  delete [] Giles_Var1;
  delete [] Giles_Var2;
  delete [] RelaxFactorAverage;
  delete [] RelaxFactorFourier;
  delete [] nSpan_iZones;

  delete [] Marker_TurboBoundIn;
  delete [] Marker_TurboBoundOut;
  delete [] Marker_Turbomachinery;
  delete [] Marker_Riemann;
  delete [] Marker_Giles;

  delete [] nBlades;
  delete [] FreeStreamTurboNormal;
}

string CConfig::GetFilename(string filename, const string& ext, int timeIter) const {

  /*--- Remove any extension --- */

  unsigned short lastindex = filename.find_last_of('.');
  filename = filename.substr(0, lastindex);

  /*--- Add the extension --- */

  filename = filename + string(ext);

  /*--- Append the zone number if multizone problems ---*/
  if (Multizone_Problem)
    filename = GetMultizone_FileName(filename, GetiZone(), ext);

  /*--- Append the zone number if multiple instance problems ---*/
  if (GetnTimeInstances() > 1)
    filename = GetMultiInstance_FileName(filename, GetiInst(), ext);

  /*--- Append the iteration number for unsteady problems ---*/
  if (GetTime_Domain()){
    filename = GetUnsteady_FileName(filename, timeIter, ext);
  }

  return filename;
}

string CConfig::GetFilename_Iter(const string& filename_iter, unsigned long curInnerIter, unsigned long curOuterIter) const {
  const auto iter = GetMultizone_Problem() ? curOuterIter : curInnerIter;
  std::stringstream iter_ss;
  iter_ss << filename_iter << "_" << std::setw(6) << std::setfill('0') << iter;
  return iter_ss.str();
}

string CConfig::GetUnsteady_FileName(string val_filename, int val_iter, const string& ext) const {

  string UnstExt, UnstFilename = std::move(val_filename);
  char buffer[50];

  /*--- Check that a positive value iteration is requested (for now). ---*/

  if (val_iter < 0) {
    SU2_MPI::Error("Requesting a negative iteration number for the restart file!!", CURRENT_FUNCTION);
  }

  unsigned short lastindex = UnstFilename.find_last_of('.');
  UnstFilename = UnstFilename.substr(0, lastindex);

  /*--- Append iteration number for unsteady cases ---*/

  if (Time_Domain) {

    if ((val_iter >= 0)    && (val_iter < 10))    SPRINTF (buffer, "_0000%d", val_iter);
    if ((val_iter >= 10)   && (val_iter < 100))   SPRINTF (buffer, "_000%d",  val_iter);
    if ((val_iter >= 100)  && (val_iter < 1000))  SPRINTF (buffer, "_00%d",   val_iter);
    if ((val_iter >= 1000) && (val_iter < 10000)) SPRINTF (buffer, "_0%d",    val_iter);
    if (val_iter >= 10000) SPRINTF (buffer, "_%d", val_iter);
    UnstExt = string(buffer);
  }
  UnstExt += ext;
  UnstFilename.append(UnstExt);

  return UnstFilename;
}

string CConfig::GetMultizone_FileName(string val_filename, int val_iZone, const string& ext) const {

    string multizone_filename = std::move(val_filename);
    char buffer[50];

    unsigned short lastindex = multizone_filename.find_last_of('.');
    multizone_filename = multizone_filename.substr(0, lastindex);

    if (Multizone_Problem) {
        SPRINTF (buffer, "_%d", SU2_TYPE::Int(val_iZone));
        multizone_filename.append(string(buffer));
    }

    multizone_filename += ext;
    return multizone_filename;
}

string CConfig::GetMultizone_HistoryFileName(string val_filename, int val_iZone, const string& ext) const {

    string multizone_filename = std::move(val_filename);
    char buffer[50];
    unsigned short lastindex = multizone_filename.find_last_of('.');
    multizone_filename = multizone_filename.substr(0, lastindex);
    if (Multizone_Problem) {
        SPRINTF (buffer, "_%d", SU2_TYPE::Int(val_iZone));
        multizone_filename.append(string(buffer));
    }
    multizone_filename += ext;
    return multizone_filename;
}

string CConfig::GetMultiInstance_FileName(string val_filename, int val_iInst, const string& ext) const {

  string multizone_filename = std::move(val_filename);
  char buffer[50];

  unsigned short lastindex = multizone_filename.find_last_of('.');
  multizone_filename = multizone_filename.substr(0, lastindex);
  SPRINTF (buffer, "_%d", SU2_TYPE::Int(val_iInst));
  multizone_filename.append(string(buffer));
  multizone_filename += ext;
  return multizone_filename;
}

string CConfig::GetMultiInstance_HistoryFileName(string val_filename, int val_iInst) const {

  string multizone_filename = std::move(val_filename);
  char buffer[50];

  unsigned short lastindex = multizone_filename.find_last_of('.');
  multizone_filename = multizone_filename.substr(0, lastindex);
  SPRINTF (buffer, "_%d", SU2_TYPE::Int(val_iInst));
  multizone_filename.append(string(buffer));

  return multizone_filename;
}

string CConfig::GetObjFunc_Extension(string val_filename) const {

  string AdjExt, Filename = std::move(val_filename);

  if (ContinuousAdjoint || DiscreteAdjoint) {

    /*--- Remove filename extension (.dat) ---*/

    unsigned short lastindex = Filename.find_last_of('.');
    Filename = Filename.substr(0, lastindex);

    if (nObj==1) {
      switch (Kind_ObjFunc[0]) {
        case DRAG_COEFFICIENT:            AdjExt = "_cd";       break;
        case LIFT_COEFFICIENT:            AdjExt = "_cl";       break;
        case SIDEFORCE_COEFFICIENT:       AdjExt = "_csf";      break;
        case INVERSE_DESIGN_PRESSURE:     AdjExt = "_invpress"; break;
        case INVERSE_DESIGN_HEATFLUX:     AdjExt = "_invheat";  break;
        case MOMENT_X_COEFFICIENT:        AdjExt = "_cmx";      break;
        case MOMENT_Y_COEFFICIENT:        AdjExt = "_cmy";      break;
        case MOMENT_Z_COEFFICIENT:        AdjExt = "_cmz";      break;
        case EFFICIENCY:                  AdjExt = "_eff";      break;
        case EQUIVALENT_AREA:             AdjExt = "_ea";       break;
        case NEARFIELD_PRESSURE:          AdjExt = "_nfp";      break;
        case FORCE_X_COEFFICIENT:         AdjExt = "_cfx";      break;
        case FORCE_Y_COEFFICIENT:         AdjExt = "_cfy";      break;
        case FORCE_Z_COEFFICIENT:         AdjExt = "_cfz";      break;
        case THRUST_COEFFICIENT:          AdjExt = "_ct";       break;
        case TORQUE_COEFFICIENT:          AdjExt = "_cq";       break;
        case TOTAL_HEATFLUX:              AdjExt = "_totheat";  break;
        case MAXIMUM_HEATFLUX:            AdjExt = "_maxheat";  break;
        case AVG_TEMPERATURE:             AdjExt = "_avtp";     break;
        case FIGURE_OF_MERIT:             AdjExt = "_merit";    break;
        case BUFFET_SENSOR:               AdjExt = "_buffet";   break;
        case SURFACE_TOTAL_PRESSURE:      AdjExt = "_pt";       break;
        case SURFACE_STATIC_PRESSURE:     AdjExt = "_pe";       break;
        case SURFACE_STATIC_TEMPERATURE:  AdjExt = "_T";        break;
        case SURFACE_MASSFLOW:            AdjExt = "_mfr";      break;
        case SURFACE_UNIFORMITY:          AdjExt = "_uniform";  break;
        case SURFACE_SECONDARY:           AdjExt = "_second";   break;
        case SURFACE_MOM_DISTORTION:      AdjExt = "_distort";  break;
        case SURFACE_SECOND_OVER_UNIFORM: AdjExt = "_sou";      break;
        case SURFACE_PRESSURE_DROP:       AdjExt = "_dp";       break;
        case SURFACE_SPECIES_0:           AdjExt = "_avgspec0"; break;
        case SURFACE_SPECIES_VARIANCE:    AdjExt = "_specvar";  break;
        case SURFACE_MACH:                AdjExt = "_mach";     break;
        case CUSTOM_OBJFUNC:              AdjExt = "_custom";   break;
        case REFERENCE_GEOMETRY:          AdjExt = "_refgeom";  break;
        case REFERENCE_NODE:              AdjExt = "_refnode";  break;
        case VOLUME_FRACTION:             AdjExt = "_volfrac";  break;
        case TOPOL_DISCRETENESS:          AdjExt = "_topdisc";  break;
        case TOPOL_COMPLIANCE:            AdjExt = "_topcomp";  break;
        case STRESS_PENALTY:              AdjExt = "_stress";   break;
      }
    }
    else{
      AdjExt = "_combo";
    }
    Filename.append(AdjExt);

    /*--- Lastly, add the .dat extension ---*/
    Filename.append(".dat");

  }

  return Filename;
}

unsigned short CConfig::GetContainerPosition(unsigned short val_eqsystem) {

  switch (val_eqsystem) {
    case RUNTIME_FLOW_SYS:      return FLOW_SOL;
    case RUNTIME_TURB_SYS:      return TURB_SOL;
    case RUNTIME_TRANS_SYS:     return TRANS_SOL;
    case RUNTIME_SPECIES_SYS:   return SPECIES_SOL;
    case RUNTIME_HEAT_SYS:      return HEAT_SOL;
    case RUNTIME_FEA_SYS:       return FEA_SOL;
    case RUNTIME_ADJFLOW_SYS:   return ADJFLOW_SOL;
    case RUNTIME_ADJTURB_SYS:   return ADJTURB_SOL;
    case RUNTIME_ADJSPECIES_SYS:return ADJSPECIES_SOL;
    case RUNTIME_ADJFEA_SYS:    return ADJFEA_SOL;
    case RUNTIME_RADIATION_SYS: return RAD_SOL;
    case RUNTIME_MULTIGRID_SYS: return 0;
  }
  return 0;
}

void CConfig::SetKind_ConvNumScheme(unsigned short val_kind_convnumscheme,
                                    CENTERED val_kind_centered, UPWIND val_kind_upwind,
                                    LIMITER val_kind_slopelimit, bool val_muscl,
                                    unsigned short val_kind_fem) {
  Kind_ConvNumScheme = val_kind_convnumscheme;
  Kind_Centered = val_kind_centered;
  Kind_Upwind = val_kind_upwind;
  Kind_FEM = val_kind_fem;
  Kind_SlopeLimit = val_kind_slopelimit;
  MUSCL = val_muscl;

}

void CConfig::SetGlobalParam(MAIN_SOLVER val_solver,
                             unsigned short val_system) {

  /*--- Set the simulation global time ---*/

  Current_UnstTime = static_cast<su2double>(TimeIter)*Delta_UnstTime;
  Current_UnstTimeND = static_cast<su2double>(TimeIter)*Delta_UnstTimeND;

  /*--- Set the solver methods ---*/

  auto SetFlowParam = [&]() {
    if (val_system == RUNTIME_FLOW_SYS) {
      SetKind_ConvNumScheme(Kind_ConvNumScheme_Flow, Kind_Centered_Flow,
                            Kind_Upwind_Flow, Kind_SlopeLimit_Flow,
                            MUSCL_Flow, NONE);
      SetKind_TimeIntScheme(Kind_TimeIntScheme_Flow);
    }
  };

  auto SetTurbParam = [&]() {
    if (val_system == RUNTIME_TURB_SYS) {
      SetKind_ConvNumScheme(Kind_ConvNumScheme_Turb, Kind_Centered_Turb,
                            Kind_Upwind_Turb, Kind_SlopeLimit_Turb,
                            MUSCL_Turb, NONE);
      SetKind_TimeIntScheme(Kind_TimeIntScheme_Turb);
    }
  };

  auto SetHeatParam = [&]() {
    if (val_system == RUNTIME_HEAT_SYS) {
      SetKind_ConvNumScheme(Kind_ConvNumScheme_Heat, Kind_Centered_Heat,
                            Kind_Upwind_Heat, Kind_SlopeLimit_Heat, MUSCL_Heat, NONE);
      SetKind_TimeIntScheme(Kind_TimeIntScheme_Heat);
    }
  };

  auto SetSpeciesParam = [&]() {
    if (val_system == RUNTIME_SPECIES_SYS) {
      SetKind_ConvNumScheme(Kind_ConvNumScheme_Species, Kind_Centered_Species,
                            Kind_Upwind_Species, Kind_SlopeLimit_Species,
                            MUSCL_Species, NONE);
      SetKind_TimeIntScheme(Kind_TimeIntScheme_Species);
    }
  };

  auto SetAdjFlowParam = [&]() {
    if (val_system == RUNTIME_ADJFLOW_SYS) {
      SetKind_ConvNumScheme(Kind_ConvNumScheme_AdjFlow, Kind_Centered_AdjFlow,
                            Kind_Upwind_AdjFlow, Kind_SlopeLimit_AdjFlow,
                            MUSCL_AdjFlow, NONE);
      SetKind_TimeIntScheme(Kind_TimeIntScheme_AdjFlow);
    }
  };

  switch (val_solver) {
    case MAIN_SOLVER::EULER: case MAIN_SOLVER::INC_EULER: case MAIN_SOLVER::NEMO_EULER:
    case MAIN_SOLVER::DISC_ADJ_EULER: case MAIN_SOLVER::DISC_ADJ_INC_EULER:
      SetFlowParam();
      break;
    case MAIN_SOLVER::NAVIER_STOKES: case MAIN_SOLVER::INC_NAVIER_STOKES: case MAIN_SOLVER::NEMO_NAVIER_STOKES:
    case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES:
      SetFlowParam();
      SetSpeciesParam();
      SetHeatParam();
      break;
    case MAIN_SOLVER::RANS: case MAIN_SOLVER::INC_RANS:
    case MAIN_SOLVER::DISC_ADJ_RANS: case MAIN_SOLVER::DISC_ADJ_INC_RANS:
      SetFlowParam();
      SetTurbParam();
      SetSpeciesParam();
      SetHeatParam();

      if (val_system == RUNTIME_TRANS_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_Turb, Kind_Centered_Turb,
                              Kind_Upwind_Turb, Kind_SlopeLimit_Turb,
                              MUSCL_Turb, NONE);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Turb);
      }
      break;
    case MAIN_SOLVER::FEM_EULER:
    case MAIN_SOLVER::FEM_NAVIER_STOKES:
    case MAIN_SOLVER::FEM_RANS:
    case MAIN_SOLVER::FEM_LES:
    case MAIN_SOLVER::DISC_ADJ_FEM_EULER:
    case MAIN_SOLVER::DISC_ADJ_FEM_NS:
    case MAIN_SOLVER::DISC_ADJ_FEM_RANS:
      if (val_system == RUNTIME_FLOW_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_FEM_Flow, Kind_Centered_Flow,
                              Kind_Upwind_Flow, Kind_SlopeLimit_Flow,
                              MUSCL_Flow, Kind_FEM_Flow);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_FEM_Flow);
      }
      break;
    case MAIN_SOLVER::ADJ_EULER:
    case MAIN_SOLVER::ADJ_NAVIER_STOKES:
      SetFlowParam();
      SetAdjFlowParam();
      break;
    case MAIN_SOLVER::ADJ_RANS:
      SetFlowParam();
      SetTurbParam();
      SetAdjFlowParam();

      if (val_system == RUNTIME_ADJTURB_SYS) {
        SetKind_ConvNumScheme(Kind_ConvNumScheme_AdjTurb, Kind_Centered_AdjTurb,
                              Kind_Upwind_AdjTurb, Kind_SlopeLimit_AdjTurb,
                              MUSCL_AdjTurb, NONE);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_AdjTurb);
      }
      break;
    case MAIN_SOLVER::HEAT_EQUATION:
    case MAIN_SOLVER::DISC_ADJ_HEAT:
      if (val_system == RUNTIME_HEAT_SYS) {
        SetKind_ConvNumScheme(NONE, CENTERED::NONE, UPWIND::NONE, LIMITER::NONE, NONE, NONE);
        SetKind_TimeIntScheme(Kind_TimeIntScheme_Heat);
      }
      break;

    case MAIN_SOLVER::FEM_ELASTICITY:
    case MAIN_SOLVER::DISC_ADJ_FEM:
      if (val_system == RUNTIME_FEA_SYS) {
        SetKind_ConvNumScheme(NONE, CENTERED::NONE, UPWIND::NONE, LIMITER::NONE, NONE, NONE);
        SetKind_TimeIntScheme(NONE);
      }
      break;

    default:
      break;
  }
}

const su2double* CConfig::GetPeriodicRotCenter(const string& val_marker) const {
  unsigned short iMarker_PerBound;
  for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++)
    if (Marker_PerBound[iMarker_PerBound] == val_marker) break;
  return Periodic_RotCenter[iMarker_PerBound];
}

const su2double* CConfig::GetPeriodicRotAngles(const string& val_marker) const {
  unsigned short iMarker_PerBound;
  for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++)
    if (Marker_PerBound[iMarker_PerBound] == val_marker) break;
  return Periodic_RotAngles[iMarker_PerBound];
}

const su2double* CConfig::GetPeriodicTranslation(const string& val_marker) const {
  unsigned short iMarker_PerBound;
  for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++)
    if (Marker_PerBound[iMarker_PerBound] == val_marker) break;
  return Periodic_Translation[iMarker_PerBound];
}

unsigned short CConfig::GetMarker_Periodic_Donor(const string& val_marker) const {
  unsigned short iMarker_PerBound, jMarker_PerBound, kMarker_All;

  /*--- Find the marker for this periodic boundary. ---*/
  for (iMarker_PerBound = 0; iMarker_PerBound < nMarker_PerBound; iMarker_PerBound++)
    if (Marker_PerBound[iMarker_PerBound] == val_marker) break;

  /*--- Find corresponding donor. ---*/
  for (jMarker_PerBound = 0; jMarker_PerBound < nMarker_PerBound; jMarker_PerBound++)
    if (Marker_PerBound[jMarker_PerBound] == Marker_PerDonor[iMarker_PerBound]) break;

  /*--- Find and return global marker index for donor boundary. ---*/
  for (kMarker_All = 0; kMarker_All < nMarker_CfgFile; kMarker_All++)
    if (Marker_PerBound[jMarker_PerBound] == Marker_All_TagBound[kMarker_All]) break;

  return kMarker_All;
}

su2double CConfig::GetActDisk_NetThrust(const string& val_marker) const {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_NetThrust[iMarker_ActDisk];
}

su2double CConfig::GetActDisk_Power(const string& val_marker) const {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_Power[iMarker_ActDisk];
}

su2double CConfig::GetActDisk_MassFlow(const string& val_marker) const {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_MassFlow[iMarker_ActDisk];
}

su2double CConfig::GetActDisk_Mach(const string& val_marker) const {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_Mach[iMarker_ActDisk];
}

su2double CConfig::GetActDisk_Force(const string& val_marker) const {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_Force[iMarker_ActDisk];
}

su2double CConfig::GetActDisk_BCThrust(const string& val_marker) const {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_BCThrust[iMarker_ActDisk];
}

su2double CConfig::GetActDisk_BCThrust_Old(const string& val_marker) const {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_BCThrust_Old[iMarker_ActDisk];
}

void CConfig::SetActDisk_BCThrust(const string& val_marker, su2double val_actdisk_bcthrust) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  ActDisk_BCThrust[iMarker_ActDisk] = val_actdisk_bcthrust;
}

void CConfig::SetActDisk_BCThrust_Old(const string& val_marker, su2double val_actdisk_bcthrust_old) {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  ActDisk_BCThrust_Old[iMarker_ActDisk] = val_actdisk_bcthrust_old;
}

su2double CConfig::GetActDisk_Area(const string& val_marker) const {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_Area[iMarker_ActDisk];
}

su2double CConfig::GetActDisk_ReverseMassFlow(const string& val_marker) const {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_ReverseMassFlow[iMarker_ActDisk];
}

su2double CConfig::GetActDisk_PressJump(const string& val_marker, unsigned short val_value) const {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_PressJump[iMarker_ActDisk][val_value];
}

su2double CConfig::GetActDisk_TempJump(const string& val_marker, unsigned short val_value) const {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_TempJump[iMarker_ActDisk][val_value];;
}

su2double CConfig::GetActDisk_Omega(const string& val_marker, unsigned short val_value) const {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if ((Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskOutlet[iMarker_ActDisk] == val_marker)) break;
  return ActDisk_Omega[iMarker_ActDisk][val_value];;
}

su2double CConfig::GetActDiskBem_CG(unsigned short iDim, string val_marker, unsigned short val_value) const {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskBemInlet_CG; iMarker_ActDisk++)
    if ((Marker_ActDiskBemInlet_CG[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskBemOutlet_CG[iMarker_ActDisk] == val_marker)) break;
  return ActDiskBem_CG[iDim][iMarker_ActDisk][val_value];
}

su2double CConfig::GetActDiskBem_Axis(unsigned short iDim, string val_marker, unsigned short val_value) const {
  unsigned short iMarker_ActDisk;
  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskBemInlet_Axis; iMarker_ActDisk++)
    if ((Marker_ActDiskBemInlet_Axis[iMarker_ActDisk] == val_marker) ||
        (Marker_ActDiskBemOutlet_Axis[iMarker_ActDisk] == val_marker)) break;
  return ActDiskBem_Axis[iDim][iMarker_ActDisk][val_value];
}

su2double CConfig::GetOutlet_MassFlow(const string& val_marker) const {
  unsigned short iMarker_Outlet;
  for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++)
    if ((Marker_Outlet[iMarker_Outlet] == val_marker)) break;
  return Outlet_MassFlow[iMarker_Outlet];
}

su2double CConfig::GetOutlet_Density(const string& val_marker) const {
  unsigned short iMarker_Outlet;
  for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++)
    if ((Marker_Outlet[iMarker_Outlet] == val_marker)) break;
  return Outlet_Density[iMarker_Outlet];
}

su2double CConfig::GetOutlet_Area(const string& val_marker) const {
  unsigned short iMarker_Outlet;
  for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++)
    if ((Marker_Outlet[iMarker_Outlet] == val_marker)) break;
  return Outlet_Area[iMarker_Outlet];
}

unsigned short CConfig::GetMarker_CfgFile_ActDiskOutlet(const string& val_marker) const {
  unsigned short iMarker_ActDisk, kMarker_All;

  /*--- Find the marker for this actuator disk inlet. ---*/

  for (iMarker_ActDisk = 0; iMarker_ActDisk < nMarker_ActDiskInlet; iMarker_ActDisk++)
    if (Marker_ActDiskInlet[iMarker_ActDisk] == val_marker) break;

  /*--- Find and return global marker index for the actuator disk outlet. ---*/

  for (kMarker_All = 0; kMarker_All < nMarker_CfgFile; kMarker_All++)
    if (Marker_ActDiskOutlet[iMarker_ActDisk] == Marker_CfgFile_TagBound[kMarker_All]) break;

  return kMarker_All;
}

unsigned short CConfig::GetMarker_CfgFile_EngineExhaust(const string& val_marker) const {
  unsigned short iMarker_Engine, kMarker_All;

  /*--- Find the marker for this engine inflow. ---*/

  for (iMarker_Engine = 0; iMarker_Engine < nMarker_EngineInflow; iMarker_Engine++)
    if (Marker_EngineInflow[iMarker_Engine] == val_marker) break;

  /*--- Find and return global marker index for the engine exhaust. ---*/

  for (kMarker_All = 0; kMarker_All < nMarker_CfgFile; kMarker_All++)
    if (Marker_EngineExhaust[iMarker_Engine] == Marker_CfgFile_TagBound[kMarker_All]) break;

  return kMarker_All;
}

bool CConfig::GetVolumetric_Movement() const {
  bool volumetric_movement = false;

  if (GetSurface_Movement(AEROELASTIC) ||
      GetSurface_Movement(AEROELASTIC_RIGID_MOTION)||
      GetSurface_Movement(EXTERNAL) ||
      GetSurface_Movement(EXTERNAL_ROTATION)){
    volumetric_movement = true;
  }

  if (Kind_SU2 == SU2_COMPONENT::SU2_DEF ||
      Kind_SU2 == SU2_COMPONENT::SU2_DOT ||
      DirectDiff)
  { volumetric_movement = true;}

  return volumetric_movement;
}

bool CConfig::GetSurface_Movement(unsigned short kind_movement) const {
  for (unsigned short iMarkerMoving = 0; iMarkerMoving < nKind_SurfaceMovement; iMarkerMoving++){
    if (Kind_SurfaceMovement[iMarkerMoving] == kind_movement){
      return true;
    }
  }
  return false;
}

unsigned short CConfig::GetMarker_Moving(const string& val_marker) const {
  unsigned short iMarker;

  /*--- Find the marker for this moving boundary. ---*/
  for (iMarker = 0; iMarker < nMarker_Moving; iMarker++)
    if (Marker_Moving[iMarker] == val_marker) break;

  return iMarker;
}

unsigned short CConfig::GetMarker_Deform_Mesh(const string& val_marker) const {
  unsigned short iMarker;

  /*--- Find the marker for this interface boundary. ---*/
  for (iMarker = 0; iMarker < nMarker_Deform_Mesh; iMarker++)
    if (Marker_Deform_Mesh[iMarker] == val_marker) break;

  return iMarker;
}

unsigned short CConfig::GetMarker_Deform_Mesh_Sym_Plane(const string& val_marker) const {
  unsigned short iMarker;

  /*--- Find the marker for this interface boundary. ---*/
  for (iMarker = 0; iMarker < nMarker_Deform_Mesh_Sym_Plane; iMarker++)
    if (Marker_Deform_Mesh_Sym_Plane[iMarker] == val_marker) break;

  return iMarker;
}

unsigned short CConfig::GetMarker_Fluid_Load(const string& val_marker) const {
  unsigned short iMarker_Fluid_Load;

  /*--- Find the marker for this interface boundary. ---*/
  for (iMarker_Fluid_Load = 0; iMarker_Fluid_Load < nMarker_Fluid_Load; iMarker_Fluid_Load++)
    if (Marker_Fluid_Load[iMarker_Fluid_Load] == val_marker) break;

  return iMarker_Fluid_Load;
}

unsigned short CConfig::GetMarker_SobolevBC(const string& val_marker) const {
  unsigned short iMarker_Sobolev;

  /*--- Find the marker for this moving boundary. ---*/
  for (iMarker_Sobolev = 0; iMarker_Sobolev < nMarker_SobolevBC; iMarker_Sobolev++)
    if (Marker_SobolevBC[iMarker_Sobolev] == val_marker) break;

  return iMarker_Sobolev;
}

su2double CConfig::GetExhaust_Temperature_Target(const string& val_marker) const {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_Temperature_Target[iMarker_EngineExhaust];
}

su2double CConfig::GetExhaust_Pressure_Target(const string& val_marker) const {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_Pressure_Target[iMarker_EngineExhaust];
}

INLET_TYPE CConfig::GetKind_Inc_Inlet(const string& val_marker) const {
  unsigned short iMarker_Inlet;
  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++)
    if (Marker_Inlet[iMarker_Inlet] == val_marker) break;
  return Kind_Inc_Inlet[iMarker_Inlet];
}

INC_OUTLET_TYPE CConfig::GetKind_Inc_Outlet(const string& val_marker) const {
  unsigned short iMarker_Outlet;
  for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++)
    if (Marker_Outlet[iMarker_Outlet] == val_marker) break;
  return Kind_Inc_Outlet[iMarker_Outlet];
}

su2double CConfig::GetInletTtotal(const string& val_marker) const {
  unsigned short iMarker_Inlet;
  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++)
    if (Marker_Inlet[iMarker_Inlet] == val_marker) break;
  return Inlet_Ttotal[iMarker_Inlet];
}

su2double CConfig::GetInletPtotal(const string& val_marker) const {
  unsigned short iMarker_Inlet;
  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++)
    if (Marker_Inlet[iMarker_Inlet] == val_marker) break;
  return Inlet_Ptotal[iMarker_Inlet];
}

void CConfig::SetInletPtotal(su2double val_pressure, const string& val_marker) {
  unsigned short iMarker_Inlet;
  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++)
    if (Marker_Inlet[iMarker_Inlet] == val_marker)
      Inlet_Ptotal[iMarker_Inlet] = val_pressure;
}

const su2double* CConfig::GetInletFlowDir(const string& val_marker) const {
  unsigned short iMarker_Inlet;
  for (iMarker_Inlet = 0; iMarker_Inlet < nMarker_Inlet; iMarker_Inlet++)
    if (Marker_Inlet[iMarker_Inlet] == val_marker) break;
  return Inlet_FlowDir[iMarker_Inlet];
}

su2double CConfig::GetInlet_Temperature(const string& val_marker) const {
  unsigned short iMarker_Supersonic_Inlet;
  for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++)
    if (Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet] == val_marker) break;
  return Inlet_Temperature[iMarker_Supersonic_Inlet];
}

su2double CConfig::GetInlet_Pressure(const string& val_marker) const {
  unsigned short iMarker_Supersonic_Inlet;
  for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++)
    if (Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet] == val_marker) break;
  return Inlet_Pressure[iMarker_Supersonic_Inlet];
}

const su2double* CConfig::GetInlet_Velocity(const string& val_marker) const {
  unsigned short iMarker_Supersonic_Inlet;
  for (iMarker_Supersonic_Inlet = 0; iMarker_Supersonic_Inlet < nMarker_Supersonic_Inlet; iMarker_Supersonic_Inlet++)
    if (Marker_Supersonic_Inlet[iMarker_Supersonic_Inlet] == val_marker) break;
  return Inlet_Velocity[iMarker_Supersonic_Inlet];
}

const su2double* CConfig::GetInlet_SpeciesVal(const string& val_marker) const {
  unsigned short iMarker_Inlet_Species;
  for (iMarker_Inlet_Species = 0; iMarker_Inlet_Species < nMarker_Inlet_Species; iMarker_Inlet_Species++)
    if (Marker_Inlet_Species[iMarker_Inlet_Species] == val_marker) break;
  return Inlet_SpeciesVal[iMarker_Inlet_Species];
}

const su2double* CConfig::GetInlet_TurbVal(const string& val_marker) const {
  /*--- If Turbulent Inlet is not provided for the marker, return free stream values. ---*/
  for (auto iMarker = 0u; iMarker < nMarker_Inlet_Turb; iMarker++) {
    if (Marker_Inlet_Turb[iMarker] == val_marker) return Inlet_TurbVal[iMarker];
  }
  if (Kind_Turb_Model == TURB_MODEL::SST) {
    return TurbIntensityAndViscRatioFreeStream;
  }
  return &NuFactor_FreeStream;
}

su2double CConfig::GetOutlet_Pressure(const string& val_marker) const {
  unsigned short iMarker_Outlet;
  for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++)
    if (Marker_Outlet[iMarker_Outlet] == val_marker) break;
  return Outlet_Pressure[iMarker_Outlet];
}

void CConfig::SetOutlet_Pressure(su2double val_pressure, const string& val_marker) {
  unsigned short iMarker_Outlet;
  for (iMarker_Outlet = 0; iMarker_Outlet < nMarker_Outlet; iMarker_Outlet++)
    if (Marker_Outlet[iMarker_Outlet] == val_marker)
      Outlet_Pressure[iMarker_Outlet] = val_pressure;
}

su2double CConfig::GetRiemann_Var1(const string& val_marker) const {
  unsigned short iMarker_Riemann;
  for (iMarker_Riemann = 0; iMarker_Riemann < nMarker_Riemann; iMarker_Riemann++)
    if (Marker_Riemann[iMarker_Riemann] == val_marker) break;
  return Riemann_Var1[iMarker_Riemann];
}

su2double CConfig::GetRiemann_Var2(const string& val_marker) const {
  unsigned short iMarker_Riemann;
  for (iMarker_Riemann = 0; iMarker_Riemann < nMarker_Riemann; iMarker_Riemann++)
    if (Marker_Riemann[iMarker_Riemann] == val_marker) break;
  return Riemann_Var2[iMarker_Riemann];
}

const su2double* CConfig::GetRiemann_FlowDir(const string& val_marker) const {
  unsigned short iMarker_Riemann;
  for (iMarker_Riemann = 0; iMarker_Riemann < nMarker_Riemann; iMarker_Riemann++)
    if (Marker_Riemann[iMarker_Riemann] == val_marker) break;
  return Riemann_FlowDir[iMarker_Riemann];
}

unsigned short CConfig::GetKind_Data_Riemann(const string& val_marker) const {
  unsigned short iMarker_Riemann;
  for (iMarker_Riemann = 0; iMarker_Riemann < nMarker_Riemann; iMarker_Riemann++)
    if (Marker_Riemann[iMarker_Riemann] == val_marker) break;
  return Kind_Data_Riemann[iMarker_Riemann];
}

su2double CConfig::GetGiles_Var1(const string& val_marker) const {
  unsigned short iMarker_Giles;
  for (iMarker_Giles = 0; iMarker_Giles < nMarker_Giles; iMarker_Giles++)
    if (Marker_Giles[iMarker_Giles] == val_marker) break;
  return Giles_Var1[iMarker_Giles];
}

void CConfig::SetGiles_Var1(su2double newVar1, const string& val_marker) {
  unsigned short iMarker_Giles;
  for (iMarker_Giles = 0; iMarker_Giles < nMarker_Giles; iMarker_Giles++)
    if (Marker_Giles[iMarker_Giles] == val_marker) break;
  Giles_Var1[iMarker_Giles] = newVar1;
}

su2double CConfig::GetGiles_Var2(const string& val_marker) const {
  unsigned short iMarker_Giles;
  for (iMarker_Giles = 0; iMarker_Giles < nMarker_Giles; iMarker_Giles++)
    if (Marker_Giles[iMarker_Giles] == val_marker) break;
  return Giles_Var2[iMarker_Giles];
}

su2double CConfig::GetGiles_RelaxFactorAverage(const string& val_marker) const {
  unsigned short iMarker_Giles;
  for (iMarker_Giles = 0; iMarker_Giles < nMarker_Giles; iMarker_Giles++)
    if (Marker_Giles[iMarker_Giles] == val_marker) break;
  return RelaxFactorAverage[iMarker_Giles];
}

su2double CConfig::GetGiles_RelaxFactorFourier(const string& val_marker) const {
  unsigned short iMarker_Giles;
  for (iMarker_Giles = 0; iMarker_Giles < nMarker_Giles; iMarker_Giles++)
    if (Marker_Giles[iMarker_Giles] == val_marker) break;
  return RelaxFactorFourier[iMarker_Giles];
}

const su2double* CConfig::GetGiles_FlowDir(const string& val_marker) const {
  unsigned short iMarker_Giles;
  for (iMarker_Giles = 0; iMarker_Giles < nMarker_Giles; iMarker_Giles++)
    if (Marker_Giles[iMarker_Giles] == val_marker) break;
  return Giles_FlowDir[iMarker_Giles];
}

unsigned short CConfig::GetKind_Data_Giles(const string& val_marker) const {
  unsigned short iMarker_Giles;
  for (iMarker_Giles = 0; iMarker_Giles < nMarker_Giles; iMarker_Giles++)
    if (Marker_Giles[iMarker_Giles] == val_marker) break;
  return Kind_Data_Giles[iMarker_Giles];
}

su2double CConfig::GetPressureOut_BC() const {
  unsigned short iMarker_BC;
  su2double pres_out = 0.0;
  for (iMarker_BC = 0; iMarker_BC < nMarker_Giles; iMarker_BC++){
    if (Kind_Data_Giles[iMarker_BC] == STATIC_PRESSURE || Kind_Data_Giles[iMarker_BC] == STATIC_PRESSURE_1D || Kind_Data_Giles[iMarker_BC] == RADIAL_EQUILIBRIUM ){
      pres_out = Giles_Var1[iMarker_BC];
    }
  }
  for (iMarker_BC = 0; iMarker_BC < nMarker_Riemann; iMarker_BC++){
    if (Kind_Data_Riemann[iMarker_BC] == STATIC_PRESSURE || Kind_Data_Riemann[iMarker_BC] == RADIAL_EQUILIBRIUM){
      pres_out = Riemann_Var1[iMarker_BC];
    }
  }
  return pres_out/Pressure_Ref;
}

void CConfig::SetPressureOut_BC(su2double val_press) {
  unsigned short iMarker_BC;
  for (iMarker_BC = 0; iMarker_BC < nMarker_Giles; iMarker_BC++){
    if (Kind_Data_Giles[iMarker_BC] == STATIC_PRESSURE || Kind_Data_Giles[iMarker_BC] == STATIC_PRESSURE_1D || Kind_Data_Giles[iMarker_BC] == RADIAL_EQUILIBRIUM ){
      Giles_Var1[iMarker_BC] = val_press*Pressure_Ref;
    }
  }
  for (iMarker_BC = 0; iMarker_BC < nMarker_Riemann; iMarker_BC++){
    if (Kind_Data_Riemann[iMarker_BC] == STATIC_PRESSURE || Kind_Data_Riemann[iMarker_BC] == RADIAL_EQUILIBRIUM){
      Riemann_Var1[iMarker_BC] = val_press*Pressure_Ref;
    }
  }
}

su2double CConfig::GetTotalPressureIn_BC() const {
  unsigned short iMarker_BC;
  su2double tot_pres_in = 0.0;
  for (iMarker_BC = 0; iMarker_BC < nMarker_Giles; iMarker_BC++){
    if (Kind_Data_Giles[iMarker_BC] == TOTAL_CONDITIONS_PT || Kind_Data_Giles[iMarker_BC] == TOTAL_CONDITIONS_PT_1D){
      tot_pres_in = Giles_Var1[iMarker_BC];
    }
  }
  for (iMarker_BC = 0; iMarker_BC < nMarker_Riemann; iMarker_BC++){
    if (Kind_Data_Riemann[iMarker_BC] == TOTAL_CONDITIONS_PT ){
      tot_pres_in = Riemann_Var1[iMarker_BC];
    }
  }
  if(nMarker_Inlet == 1 && Kind_Inlet == INLET_TYPE::TOTAL_CONDITIONS){
    tot_pres_in = Inlet_Ptotal[0];
  }
  return tot_pres_in/Pressure_Ref;
}

su2double CConfig::GetTotalTemperatureIn_BC() const {
  unsigned short iMarker_BC;
  su2double tot_temp_in = 0.0;
  for (iMarker_BC = 0; iMarker_BC < nMarker_Giles; iMarker_BC++){
    if (Kind_Data_Giles[iMarker_BC] == TOTAL_CONDITIONS_PT || Kind_Data_Giles[iMarker_BC] == TOTAL_CONDITIONS_PT_1D){
      tot_temp_in = Giles_Var2[iMarker_BC];
    }
  }
  for (iMarker_BC = 0; iMarker_BC < nMarker_Riemann; iMarker_BC++){
    if (Kind_Data_Riemann[iMarker_BC] == TOTAL_CONDITIONS_PT ){
      tot_temp_in = Riemann_Var2[iMarker_BC];
    }
  }

  if(nMarker_Inlet == 1 && Kind_Inlet == INLET_TYPE::TOTAL_CONDITIONS){
    tot_temp_in = Inlet_Ttotal[0];
  }
  return tot_temp_in/Temperature_Ref;
}

void CConfig::SetTotalTemperatureIn_BC(su2double val_temp) {
  unsigned short iMarker_BC;
  for (iMarker_BC = 0; iMarker_BC < nMarker_Giles; iMarker_BC++){
    if (Kind_Data_Giles[iMarker_BC] == TOTAL_CONDITIONS_PT || Kind_Data_Giles[iMarker_BC] == TOTAL_CONDITIONS_PT_1D){
      Giles_Var2[iMarker_BC] = val_temp*Temperature_Ref;
    }
  }
  for (iMarker_BC = 0; iMarker_BC < nMarker_Riemann; iMarker_BC++){
    if (Kind_Data_Riemann[iMarker_BC] == TOTAL_CONDITIONS_PT ){
      Riemann_Var2[iMarker_BC] = val_temp*Temperature_Ref;
    }
  }

  if(nMarker_Inlet == 1 && Kind_Inlet == INLET_TYPE::TOTAL_CONDITIONS){
    Inlet_Ttotal[0] = val_temp*Temperature_Ref;
  }
}

su2double CConfig::GetFlowAngleIn_BC() const {
  unsigned short iMarker_BC;
  su2double alpha_in = 0.0;
  for (iMarker_BC = 0; iMarker_BC < nMarker_Giles; iMarker_BC++){
    if (Kind_Data_Giles[iMarker_BC] == TOTAL_CONDITIONS_PT || Kind_Data_Giles[iMarker_BC] == TOTAL_CONDITIONS_PT_1D){
      alpha_in = atan(Giles_FlowDir[iMarker_BC][1]/Giles_FlowDir[iMarker_BC][0]);
    }
  }
  for (iMarker_BC = 0; iMarker_BC < nMarker_Riemann; iMarker_BC++){
    if (Kind_Data_Riemann[iMarker_BC] == TOTAL_CONDITIONS_PT ){
      alpha_in = atan(Riemann_FlowDir[iMarker_BC][1]/Riemann_FlowDir[iMarker_BC][0]);
    }
  }

  if(nMarker_Inlet == 1 && Kind_Inlet == INLET_TYPE::TOTAL_CONDITIONS){
    alpha_in = atan(Inlet_FlowDir[0][1]/Inlet_FlowDir[0][0]);
  }

  return alpha_in;
}

su2double CConfig::GetIncInlet_BC() const {

  su2double val_out = 0.0;

  if (nMarker_Inlet > 0) {
    if (Kind_Inc_Inlet[0] == INLET_TYPE::VELOCITY_INLET)
      val_out = Inlet_Ptotal[0]/Velocity_Ref;
    else if (Kind_Inc_Inlet[0] == INLET_TYPE::PRESSURE_INLET)
      val_out = Inlet_Ptotal[0]/Pressure_Ref;
  }

  return val_out;
}

void CConfig::SetIncInlet_BC(su2double val_in) {

  if (nMarker_Inlet > 0) {
    if (Kind_Inc_Inlet[0] == INLET_TYPE::VELOCITY_INLET)
      Inlet_Ptotal[0] = val_in*Velocity_Ref;
    else if (Kind_Inc_Inlet[0] == INLET_TYPE::PRESSURE_INLET)
      Inlet_Ptotal[0] = val_in*Pressure_Ref;
  }
}

su2double CConfig::GetIncTemperature_BC() const {

  su2double val_out = 0.0;

  if (nMarker_Inlet > 0)
    val_out = Inlet_Ttotal[0]/Temperature_Ref;

  return val_out;
}

void CConfig::SetIncTemperature_BC(su2double val_temperature) {
  if (nMarker_Inlet > 0)
    Inlet_Ttotal[0] = val_temperature*Temperature_Ref;
}

su2double CConfig::GetIncPressureOut_BC() const {

  su2double pressure_out = 0.0;

  if (nMarker_FarField > 0){
    pressure_out = Pressure_FreeStreamND;
  } else if (nMarker_Outlet > 0) {
    pressure_out = Outlet_Pressure[0]/Pressure_Ref;
  }

  return pressure_out;
}

void CConfig::SetIncPressureOut_BC(su2double val_pressure) {

  if (nMarker_FarField > 0){
    Pressure_FreeStreamND = val_pressure;
  } else if (nMarker_Outlet > 0) {
    Outlet_Pressure[0] = val_pressure*Pressure_Ref;
  }

}

su2double CConfig::GetIsothermal_Temperature(const string& val_marker) const {

  for (unsigned short iMarker_Isothermal = 0; iMarker_Isothermal < nMarker_Isothermal; iMarker_Isothermal++)
    if (Marker_Isothermal[iMarker_Isothermal] == val_marker)
      return Isothermal_Temperature[iMarker_Isothermal];

  return Isothermal_Temperature[0];
}

su2double CConfig::GetWall_HeatFlux(const string& val_marker) const {

  for (unsigned short iMarker_HeatFlux = 0; iMarker_HeatFlux < nMarker_HeatFlux; iMarker_HeatFlux++)
    if (Marker_HeatFlux[iMarker_HeatFlux] == val_marker)
      return Heat_Flux[iMarker_HeatFlux];

  return Heat_Flux[0];
}

su2double CConfig::GetWall_HeatTransfer_Coefficient(const string& val_marker) const {

  for (unsigned short iMarker_HeatTransfer = 0; iMarker_HeatTransfer < nMarker_HeatTransfer; iMarker_HeatTransfer++)
    if (Marker_HeatTransfer[iMarker_HeatTransfer] == val_marker)
      return HeatTransfer_Coeff[iMarker_HeatTransfer];

  return HeatTransfer_Coeff[0];
}

su2double CConfig::GetWall_HeatTransfer_Temperature(const string& val_marker) const {

  for (unsigned short iMarker_HeatTransfer = 0; iMarker_HeatTransfer < nMarker_HeatTransfer; iMarker_HeatTransfer++)
    if (Marker_HeatTransfer[iMarker_HeatTransfer] == val_marker)
      return HeatTransfer_WallTemp[iMarker_HeatTransfer];

  return HeatTransfer_WallTemp[0];
}

pair<WALL_TYPE, su2double> CConfig::GetWallRoughnessProperties(const string& val_marker) const {
  su2double roughness = 0.0;
  for (auto iMarker = 0u; iMarker < nRough_Wall; iMarker++) {
    if (val_marker.compare(Marker_RoughWall[iMarker]) == 0) {
      roughness = Roughness_Height[iMarker];
      break;
    }
  }
  return make_pair(roughness > 0 ? WALL_TYPE::ROUGH : WALL_TYPE::SMOOTH, roughness);
}

WALL_FUNCTIONS CConfig::GetWallFunction_Treatment(const string& val_marker) const {

  WALL_FUNCTIONS WallFunction = WALL_FUNCTIONS::NONE;

  for(unsigned short iMarker=0; iMarker<nMarker_WallFunctions; iMarker++) {
    if(Marker_WallFunctions[iMarker] == val_marker) {
      WallFunction = Kind_WallFunctions[iMarker];
      break;
    }
  }

  return WallFunction;
}

const unsigned short* CConfig::GetWallFunction_IntInfo(const string& val_marker) const {
  unsigned short *intInfo = nullptr;

  for(unsigned short iMarker=0; iMarker<nMarker_WallFunctions; iMarker++) {
    if(Marker_WallFunctions[iMarker] == val_marker) {
      intInfo = IntInfo_WallFunctions[iMarker];
      break;
    }
  }

  return intInfo;
}

const su2double* CConfig::GetWallFunction_DoubleInfo(const string& val_marker) const {
  su2double *doubleInfo = nullptr;

  for(unsigned short iMarker=0; iMarker<nMarker_WallFunctions; iMarker++) {
    if(Marker_WallFunctions[iMarker] == val_marker) {
      doubleInfo = DoubleInfo_WallFunctions[iMarker];
      break;
    }
  }

  return doubleInfo;
}

su2double CConfig::GetEngineInflow_Target(const string& val_marker) const {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return EngineInflow_Target[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_Pressure(const string& val_marker) const {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_Pressure[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_MassFlow(const string& val_marker) const {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_MassFlow[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_ReverseMassFlow(const string& val_marker) const {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_ReverseMassFlow[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_TotalPressure(const string& val_marker) const {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_TotalPressure[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_Temperature(const string& val_marker) const {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_Temperature[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_TotalTemperature(const string& val_marker) const {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_TotalTemperature[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_RamDrag(const string& val_marker) const {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_RamDrag[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_Force(const string& val_marker) const {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_Force[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_Power(const string& val_marker) const {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_Power[iMarker_EngineInflow];
}

su2double CConfig::GetInflow_Mach(const string& val_marker) const {
  unsigned short iMarker_EngineInflow;
  for (iMarker_EngineInflow = 0; iMarker_EngineInflow < nMarker_EngineInflow; iMarker_EngineInflow++)
    if (Marker_EngineInflow[iMarker_EngineInflow] == val_marker) break;
  return Inflow_Mach[iMarker_EngineInflow];
}

su2double CConfig::GetExhaust_Pressure(const string& val_marker) const {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_Pressure[iMarker_EngineExhaust];
}

su2double CConfig::GetExhaust_Temperature(const string& val_marker) const {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_Temperature[iMarker_EngineExhaust];
}

su2double CConfig::GetExhaust_MassFlow(const string& val_marker) const {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_MassFlow[iMarker_EngineExhaust];
}

su2double CConfig::GetExhaust_TotalPressure(const string& val_marker) const {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_TotalPressure[iMarker_EngineExhaust];
}

su2double CConfig::GetExhaust_TotalTemperature(const string& val_marker) const {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_TotalTemperature[iMarker_EngineExhaust];
}

su2double CConfig::GetExhaust_GrossThrust(const string& val_marker) const {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_GrossThrust[iMarker_EngineExhaust];
}

su2double CConfig::GetExhaust_Force(const string& val_marker) const {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_Force[iMarker_EngineExhaust];
}

su2double CConfig::GetExhaust_Power(const string& val_marker) const {
  unsigned short iMarker_EngineExhaust;
  for (iMarker_EngineExhaust = 0; iMarker_EngineExhaust < nMarker_EngineExhaust; iMarker_EngineExhaust++)
    if (Marker_EngineExhaust[iMarker_EngineExhaust] == val_marker) break;
  return Exhaust_Power[iMarker_EngineExhaust];
}

su2double CConfig::GetActDiskInlet_Pressure(const string& val_marker) const {
  unsigned short iMarker_ActDiskInlet;
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++)
    if (Marker_ActDiskInlet[iMarker_ActDiskInlet] == val_marker) break;
  return ActDiskInlet_Pressure[iMarker_ActDiskInlet];
}

su2double CConfig::GetActDiskInlet_TotalPressure(const string& val_marker) const {
  unsigned short iMarker_ActDiskInlet;
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++)
    if (Marker_ActDiskInlet[iMarker_ActDiskInlet] == val_marker) break;
  return ActDiskInlet_TotalPressure[iMarker_ActDiskInlet];
}

su2double CConfig::GetActDiskInlet_RamDrag(const string& val_marker) const {
  unsigned short iMarker_ActDiskInlet;
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++)
    if (Marker_ActDiskInlet[iMarker_ActDiskInlet] == val_marker) break;
  return ActDiskInlet_RamDrag[iMarker_ActDiskInlet];
}

su2double CConfig::GetActDiskInlet_Force(const string& val_marker) const {
  unsigned short iMarker_ActDiskInlet;
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++)
    if (Marker_ActDiskInlet[iMarker_ActDiskInlet] == val_marker) break;
  return ActDiskInlet_Force[iMarker_ActDiskInlet];
}

su2double CConfig::GetActDiskInlet_Power(const string& val_marker) const {
  unsigned short iMarker_ActDiskInlet;
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++)
    if (Marker_ActDiskInlet[iMarker_ActDiskInlet] == val_marker) break;
  return ActDiskInlet_Power[iMarker_ActDiskInlet];
}

su2double CConfig::GetActDiskOutlet_Pressure(const string& val_marker) const {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_Pressure[iMarker_ActDiskOutlet];
}

su2double CConfig::GetActDiskOutlet_TotalPressure(const string& val_marker) const {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_TotalPressure[iMarker_ActDiskOutlet];
}

su2double CConfig::GetActDiskOutlet_GrossThrust(const string& val_marker) const {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_GrossThrust[iMarker_ActDiskOutlet];
}

su2double CConfig::GetActDiskOutlet_Force(const string& val_marker) const {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_Force[iMarker_ActDiskOutlet];
}

su2double CConfig::GetActDiskOutlet_Power(const string& val_marker) const {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_Power[iMarker_ActDiskOutlet];
}

su2double CConfig::GetActDiskOutlet_Thrust_BEM(string val_marker) const {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_Thrust_BEM[iMarker_ActDiskOutlet];
}

su2double CConfig::GetActDiskOutlet_Torque_BEM(string val_marker) const {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_Torque_BEM[iMarker_ActDiskOutlet];
}

su2double CConfig::GetActDiskInlet_Temperature(const string& val_marker) const {
  unsigned short iMarker_ActDiskInlet;
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++)
    if (Marker_ActDiskInlet[iMarker_ActDiskInlet] == val_marker) break;
  return ActDiskInlet_Temperature[iMarker_ActDiskInlet];
}

su2double CConfig::GetActDiskInlet_TotalTemperature(const string& val_marker) const {
  unsigned short iMarker_ActDiskInlet;
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++)
    if (Marker_ActDiskInlet[iMarker_ActDiskInlet] == val_marker) break;
  return ActDiskInlet_TotalTemperature[iMarker_ActDiskInlet];
}

su2double CConfig::GetActDiskOutlet_Temperature(const string& val_marker) const {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_Temperature[iMarker_ActDiskOutlet];
}

su2double CConfig::GetActDiskOutlet_TotalTemperature(const string& val_marker) const {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_TotalTemperature[iMarker_ActDiskOutlet];
}

su2double CConfig::GetActDiskInlet_MassFlow(const string& val_marker) const {
  unsigned short iMarker_ActDiskInlet;
  for (iMarker_ActDiskInlet = 0; iMarker_ActDiskInlet < nMarker_ActDiskInlet; iMarker_ActDiskInlet++)
    if (Marker_ActDiskInlet[iMarker_ActDiskInlet] == val_marker) break;
  return ActDiskInlet_MassFlow[iMarker_ActDiskInlet];
}

su2double CConfig::GetActDiskOutlet_MassFlow(const string& val_marker) const {
  unsigned short iMarker_ActDiskOutlet;
  for (iMarker_ActDiskOutlet = 0; iMarker_ActDiskOutlet < nMarker_ActDiskOutlet; iMarker_ActDiskOutlet++)
    if (Marker_ActDiskOutlet[iMarker_ActDiskOutlet] == val_marker) break;
  return ActDiskOutlet_MassFlow[iMarker_ActDiskOutlet];
}

su2double CConfig::GetDispl_Value(const string& val_marker) const {
  unsigned short iMarker_Displacement;
  for (iMarker_Displacement = 0; iMarker_Displacement < nMarker_Displacement; iMarker_Displacement++)
    if (Marker_Displacement[iMarker_Displacement] == val_marker) break;
  return Displ_Value[iMarker_Displacement];
}

su2double CConfig::GetLoad_Value(const string& val_marker) const {
  unsigned short iMarker_Load;
  for (iMarker_Load = 0; iMarker_Load < nMarker_Load; iMarker_Load++)
    if (Marker_Load[iMarker_Load] == val_marker) break;
  return Load_Value[iMarker_Load];
}

su2double CConfig::GetDamper_Constant(const string& val_marker) const {
  unsigned short iMarker_Damper;
  for (iMarker_Damper = 0; iMarker_Damper < nMarker_Damper; iMarker_Damper++)
    if (Marker_Damper[iMarker_Damper] == val_marker) break;
  return Damper_Constant[iMarker_Damper];
}

su2double CConfig::GetLoad_Dir_Value(const string& val_marker) const {
  unsigned short iMarker_Load_Dir;
  for (iMarker_Load_Dir = 0; iMarker_Load_Dir < nMarker_Load_Dir; iMarker_Load_Dir++)
    if (Marker_Load_Dir[iMarker_Load_Dir] == val_marker) break;
  return Load_Dir_Value[iMarker_Load_Dir];
}

su2double CConfig::GetLoad_Dir_Multiplier(const string& val_marker) const {
  unsigned short iMarker_Load_Dir;
  for (iMarker_Load_Dir = 0; iMarker_Load_Dir < nMarker_Load_Dir; iMarker_Load_Dir++)
    if (Marker_Load_Dir[iMarker_Load_Dir] == val_marker) break;
  return Load_Dir_Multiplier[iMarker_Load_Dir];
}

su2double CConfig::GetDisp_Dir_Value(const string& val_marker) const {
  unsigned short iMarker_Disp_Dir;
  for (iMarker_Disp_Dir = 0; iMarker_Disp_Dir < nMarker_Disp_Dir; iMarker_Disp_Dir++)
    if (Marker_Disp_Dir[iMarker_Disp_Dir] == val_marker) break;
  return Disp_Dir_Value[iMarker_Disp_Dir];
}

su2double CConfig::GetDisp_Dir_Multiplier(const string& val_marker) const {
  unsigned short iMarker_Disp_Dir;
  for (iMarker_Disp_Dir = 0; iMarker_Disp_Dir < nMarker_Disp_Dir; iMarker_Disp_Dir++)
    if (Marker_Disp_Dir[iMarker_Disp_Dir] == val_marker) break;
  return Disp_Dir_Multiplier[iMarker_Disp_Dir];
}

const su2double* CConfig::GetLoad_Dir(const string& val_marker) const {
  unsigned short iMarker_Load_Dir;
  for (iMarker_Load_Dir = 0; iMarker_Load_Dir < nMarker_Load_Dir; iMarker_Load_Dir++)
    if (Marker_Load_Dir[iMarker_Load_Dir] == val_marker) break;
  return Load_Dir[iMarker_Load_Dir];
}

const su2double* CConfig::GetDisp_Dir(const string& val_marker) const {
  unsigned short iMarker_Disp_Dir;
  for (iMarker_Disp_Dir = 0; iMarker_Disp_Dir < nMarker_Disp_Dir; iMarker_Disp_Dir++)
    if (Marker_Disp_Dir[iMarker_Disp_Dir] == val_marker) break;
  return Disp_Dir[iMarker_Disp_Dir];
}

su2double CConfig::GetWall_Emissivity(const string& val_marker) const {
  for (auto iMarker = 0u; iMarker < nMarker_Emissivity; iMarker++)
    if (Marker_Emissivity[iMarker] == val_marker)
      return Wall_Emissivity[iMarker];
  return 0;
}

bool CConfig::GetMarker_StrongBC(const string& val_marker) const {

  unsigned short iMarker_StrongBC = 0;

  for (iMarker_StrongBC = 0; iMarker_StrongBC < nMarker_StrongBC; iMarker_StrongBC++)
    if (Marker_StrongBC[iMarker_StrongBC] == val_marker) return true;

  return false;
}

short CConfig::FindInterfaceMarker(unsigned short iInterface) const {

  /*--- The names of the two markers that form the interface. ---*/
  const auto& sideA = Marker_ZoneInterface[2*iInterface];
  const auto& sideB = Marker_ZoneInterface[2*iInterface+1];
  for (unsigned short iMarker = 0; iMarker < nMarker_All; iMarker++) {
    /*--- If the marker is sideA or sideB of the interface (order does not matter). ---*/
    const auto& tag = Marker_All_TagBound[iMarker];
    if ((tag == sideA) || (tag == sideB)) return iMarker;
  }
  return -1;
}

void CConfig::Tick(double *val_start_time) {

#ifdef PROFILE
  *val_start_time = SU2_MPI::Wtime();
#endif

}

void CConfig::Tock(double val_start_time, const string& val_function_name, int val_group_id) {

#ifdef PROFILE

  double val_stop_time = 0.0, val_elapsed_time = 0.0;

  val_stop_time = SU2_MPI::Wtime();

  /*--- Compute the elapsed time for this subroutine ---*/
  val_elapsed_time = val_stop_time - val_start_time;

  /*--- Store the subroutine name and the elapsed time ---*/
  Profile_Function_tp.push_back(val_function_name);
  Profile_Time_tp.push_back(val_elapsed_time);
  Profile_ID_tp.push_back(val_group_id);

#endif

}

void CConfig::SetProfilingCSV() {

#ifdef PROFILE

  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
#ifdef HAVE_MPI
  SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);
  SU2_MPI::Comm_size(SU2_MPI::GetComm(), &size);
#endif

  /*--- Each rank has the same stack trace, so the they have the same
   function calls and ordering in the vectors. We're going to reduce
   the timings from each rank and extract the avg, min, and max timings. ---*/

  /*--- First, create a local mapping, so that we can extract the
   min and max values for each function. ---*/

  for (unsigned int i = 0; i < Profile_Function_tp.size(); i++) {

    /*--- Add the function and initialize if not already stored (the ID
     only needs to be stored the first time).---*/
    if (Profile_Map_tp.find(Profile_Function_tp[i]) == Profile_Map_tp.end()) {

      vector<int> profile; profile.push_back(i);
      Profile_Map_tp.insert(pair<string,vector<int> >(Profile_Function_tp[i],profile));

    } else {

      /*--- This function has already been added, so simply increment the
       number of calls and total time for this function. ---*/

      Profile_Map_tp[Profile_Function_tp[i]].push_back(i);

    }
  }

  /*--- We now have everything gathered by function name, so we can loop over
   each function and store the min/max times. ---*/

  int map_size = 0;
  for (map<string,vector<int> >::iterator it=Profile_Map_tp.begin(); it!=Profile_Map_tp.end(); ++it) {
    map_size++;
  }

  /*--- Allocate and initialize memory ---*/

  double *l_min_red = NULL, *l_max_red = NULL, *l_tot_red = NULL, *l_avg_red = NULL;
  int *n_calls_red = NULL;
  double* l_min = new double[map_size];
  double* l_max = new double[map_size];
  double* l_tot = new double[map_size];
  double* l_avg = new double[map_size];
  int* n_calls  = new int[map_size];
  for (int i = 0; i < map_size; i++)
  {
    l_min[i]   = 1e10;
    l_max[i]   = 0.0;
    l_tot[i]   = 0.0;
    l_avg[i]   = 0.0;
    n_calls[i] = 0;
  }

  /*--- Collect the info for each function from the current rank ---*/

  int func_counter = 0;
  for (map<string,vector<int> >::iterator it=Profile_Map_tp.begin(); it!=Profile_Map_tp.end(); ++it) {

    for (unsigned int i = 0; i < (it->second).size(); i++) {
      n_calls[func_counter]++;
      l_tot[func_counter] += Profile_Time_tp[(it->second)[i]];
      if (Profile_Time_tp[(it->second)[i]] < l_min[func_counter])
        l_min[func_counter] = Profile_Time_tp[(it->second)[i]];
      if (Profile_Time_tp[(it->second)[i]] > l_max[func_counter])
        l_max[func_counter] = Profile_Time_tp[(it->second)[i]];

    }
    l_avg[func_counter] = l_tot[func_counter]/((double)n_calls[func_counter]);
    func_counter++;
  }

  /*--- Now reduce the data ---*/

  if (rank == MASTER_NODE) {
    l_min_red = new double[map_size];
    l_max_red = new double[map_size];
    l_tot_red = new double[map_size];
    l_avg_red = new double[map_size];
    n_calls_red  = new int[map_size];
  }

#ifdef HAVE_MPI
  MPI_Reduce(n_calls, n_calls_red, map_size, MPI_INT, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());
  MPI_Reduce(l_tot, l_tot_red, map_size, MPI_DOUBLE, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());
  MPI_Reduce(l_avg, l_avg_red, map_size, MPI_DOUBLE, MPI_SUM, MASTER_NODE, SU2_MPI::GetComm());
  MPI_Reduce(l_min, l_min_red, map_size, MPI_DOUBLE, MPI_MIN, MASTER_NODE, SU2_MPI::GetComm());
  MPI_Reduce(l_max, l_max_red, map_size, MPI_DOUBLE, MPI_MAX, MASTER_NODE, SU2_MPI::GetComm());
#else
  memcpy(n_calls_red, n_calls, map_size*sizeof(int));
  memcpy(l_tot_red,   l_tot,   map_size*sizeof(double));
  memcpy(l_avg_red,   l_avg,   map_size*sizeof(double));
  memcpy(l_min_red,   l_min,   map_size*sizeof(double));
  memcpy(l_max_red,   l_max,   map_size*sizeof(double));
#endif

  /*--- The master rank will write the file ---*/

  if (rank == MASTER_NODE) {

    /*--- Take averages over all ranks on the master ---*/

    for (int i = 0; i < map_size; i++) {
      l_tot_red[i]   = l_tot_red[i]/(double)size;
      l_avg_red[i]   = l_avg_red[i]/(double)size;
      n_calls_red[i] = n_calls_red[i]/size;
    }

    /*--- Now write a CSV file with the processed results ---*/

    ofstream Profile_File;
    Profile_File.precision(15);
    Profile_File.open("profiling.csv");

    /*--- Create the CSV header ---*/

    Profile_File << "\"Function_Name\", \"N_Calls\", \"Avg_Total_Time\", \"Avg_Time\", \"Min_Time\", \"Max_Time\", \"Function_ID\"" << endl;

    /*--- Loop through the map and write the results to the file ---*/

    func_counter = 0;
    for (map<string,vector<int> >::iterator it=Profile_Map_tp.begin(); it!=Profile_Map_tp.end(); ++it) {

      Profile_File << scientific << it->first << ", " << n_calls_red[func_counter] << ", " << l_tot_red[func_counter] << ", " << l_avg_red[func_counter] << ", " << l_min_red[func_counter] << ", " << l_max_red[func_counter] << ", " << (int)Profile_ID_tp[(it->second)[0]] << endl;
      func_counter++;
    }

    Profile_File.close();

  }

  delete [] l_min;
  delete [] l_max;
  delete [] l_avg;
  delete [] l_tot;
  delete [] n_calls;
  if (rank == MASTER_NODE) {
    delete [] l_min_red;
    delete [] l_max_red;
    delete [] l_avg_red;
    delete [] l_tot_red;
    delete [] n_calls_red;
  }

#endif

}

void CConfig::GEMM_Tick(double *val_start_time) const {

#ifdef PROFILE

#ifdef HAVE_MKL
  *val_start_time = dsecnd();
#else
  *val_start_time = SU2_MPI::Wtime();
#endif

#endif

}

void CConfig::GEMM_Tock(double val_start_time, int M, int N, int K) const {

#ifdef PROFILE

  /* Determine the timing value. The actual function called depends on
     the type of executable. */
  double val_stop_time = 0.0;

#ifdef HAVE_MKL
  val_stop_time = dsecnd();
#else
  val_stop_time = SU2_MPI::Wtime();
#endif

  /* Compute the elapsed time. */
  const double val_elapsed_time = val_stop_time - val_start_time;

  /* Create the CLong3T from the M-N-K values and check if it is already
     stored in the map GEMM_Profile_MNK. */
  CLong3T MNK(M, N, K);
  map<CLong3T, int>::iterator MI = GEMM_Profile_MNK.find(MNK);

  if(MI == GEMM_Profile_MNK.end()) {

    /* Entry is not present yet. Create it. */
    const int ind = GEMM_Profile_MNK.size();
    GEMM_Profile_MNK[MNK] = ind;

    GEMM_Profile_NCalls.push_back(1);
    GEMM_Profile_TotTime.push_back(val_elapsed_time);
    GEMM_Profile_MinTime.push_back(val_elapsed_time);
    GEMM_Profile_MaxTime.push_back(val_elapsed_time);
  }
  else {

    /* Entry is already present. Determine its index in the
       map and update the corresponding vectors. */
    const int ind = MI->second;
    ++GEMM_Profile_NCalls[ind];
    GEMM_Profile_TotTime[ind] += val_elapsed_time;
    GEMM_Profile_MinTime[ind]  = min(GEMM_Profile_MinTime[ind], val_elapsed_time);
    GEMM_Profile_MaxTime[ind]  = max(GEMM_Profile_MaxTime[ind], val_elapsed_time);
  }

#endif

}

void CConfig::GEMMProfilingCSV() {

#ifdef PROFILE

  /* Initialize the rank to the master node. */
  int rank = MASTER_NODE;

#ifdef HAVE_MPI
  /* Parallel executable. The profiling data must be sent to the master node.
     First determine the rank and size. */
  int size;
  SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);
  SU2_MPI::Comm_size(SU2_MPI::GetComm(), &size);

  /* Check for the master node. */
  if(rank == MASTER_NODE) {

    /* Master node. Loop over the other ranks to receive their data. */
    for(int proc=1; proc<size; ++proc) {

      /* Block until a message from this processor arrives. Determine
         the number of entries in the receive buffers. */
      SU2_MPI::Status status;
      SU2_MPI::Probe(proc, 0, SU2_MPI::GetComm(), &status);

      int nEntries;
      SU2_MPI::Get_count(&status, MPI_LONG, &nEntries);

      /* Allocate the memory for the receive buffers and receive the
         three messages using blocking receives. */
      vector<long>   recvBufNCalls(nEntries);
      vector<double> recvBufTotTime(nEntries);
      vector<double> recvBufMinTime(nEntries);
      vector<double> recvBufMaxTime(nEntries);
      vector<long>   recvBufMNK(3*nEntries);

      SU2_MPI::Recv(recvBufNCalls.data(), recvBufNCalls.size(),
                    MPI_LONG, proc, 0, SU2_MPI::GetComm(), &status);
      SU2_MPI::Recv(recvBufTotTime.data(), recvBufTotTime.size(),
                    MPI_DOUBLE, proc, 1, SU2_MPI::GetComm(), &status);
      SU2_MPI::Recv(recvBufMinTime.data(), recvBufMinTime.size(),
                    MPI_DOUBLE, proc, 2, SU2_MPI::GetComm(), &status);
      SU2_MPI::Recv(recvBufMaxTime.data(), recvBufMaxTime.size(),
                    MPI_DOUBLE, proc, 3, SU2_MPI::GetComm(), &status);
      SU2_MPI::Recv(recvBufMNK.data(), recvBufMNK.size(),
                    MPI_LONG, proc, 4, SU2_MPI::GetComm(), &status);

      /* Loop over the number of entries. */
      for(int i=0; i<nEntries; ++i) {

        /* Create the CLong3T from the M-N-K values and check if it is already
           stored in the map GEMM_Profile_MNK. */
        CLong3T MNK(recvBufMNK[3*i], recvBufMNK[3*i+1], recvBufMNK[3*i+2]);
        map<CLong3T, int>::iterator MI = GEMM_Profile_MNK.find(MNK);

        if(MI == GEMM_Profile_MNK.end()) {

          /* Entry is not present yet. Create it. */
          const int ind = GEMM_Profile_MNK.size();
          GEMM_Profile_MNK[MNK] = ind;

          GEMM_Profile_NCalls.push_back(recvBufNCalls[i]);
          GEMM_Profile_TotTime.push_back(recvBufTotTime[i]);
          GEMM_Profile_MinTime.push_back(recvBufMinTime[i]);
          GEMM_Profile_MaxTime.push_back(recvBufMaxTime[i]);
        }
        else {

          /* Entry is already present. Determine its index in the
             map and update the corresponding vectors. */
          const int ind = MI->second;
          GEMM_Profile_NCalls[ind]  += recvBufNCalls[i];
          GEMM_Profile_TotTime[ind] += recvBufTotTime[i];
          GEMM_Profile_MinTime[ind]  = min(GEMM_Profile_MinTime[ind], recvBufMinTime[i]);
          GEMM_Profile_MaxTime[ind]  = max(GEMM_Profile_MaxTime[ind], recvBufMaxTime[i]);
        }
      }
    }
  }
  else {

    /* Not the master node. Create the send buffer for the MNK data. */
    vector<long> sendBufMNK(3*GEMM_Profile_NCalls.size());
    for(map<CLong3T, int>::iterator MI =GEMM_Profile_MNK.begin();
                                    MI!=GEMM_Profile_MNK.end(); ++MI) {

      const int ind = 3*MI->second;
      sendBufMNK[ind]   = MI->first.long0;
      sendBufMNK[ind+1] = MI->first.long1;
      sendBufMNK[ind+2] = MI->first.long2;
    }

    /* Send the data to the master node using blocking sends. */
    SU2_MPI::Send(GEMM_Profile_NCalls.data(), GEMM_Profile_NCalls.size(),
                  MPI_LONG, MASTER_NODE, 0, SU2_MPI::GetComm());
    SU2_MPI::Send(GEMM_Profile_TotTime.data(), GEMM_Profile_TotTime.size(),
                  MPI_DOUBLE, MASTER_NODE, 1, SU2_MPI::GetComm());
    SU2_MPI::Send(GEMM_Profile_MinTime.data(), GEMM_Profile_MinTime.size(),
                  MPI_DOUBLE, MASTER_NODE, 2, SU2_MPI::GetComm());
    SU2_MPI::Send(GEMM_Profile_MaxTime.data(), GEMM_Profile_MaxTime.size(),
                  MPI_DOUBLE, MASTER_NODE, 3, SU2_MPI::GetComm());
    SU2_MPI::Send(sendBufMNK.data(), sendBufMNK.size(),
                  MPI_LONG, MASTER_NODE, 4, SU2_MPI::GetComm());
  }

#endif

  /*--- The master rank will write the file ---*/
  if (rank == MASTER_NODE) {

    /* Store the elements of the map GEMM_Profile_MNK in
       vectors for post processing reasons. */
    const unsigned int nItems = GEMM_Profile_MNK.size();
    vector<long> M(nItems), N(nItems), K(nItems);
    for(map<CLong3T, int>::iterator MI =GEMM_Profile_MNK.begin();
                                    MI!=GEMM_Profile_MNK.end(); ++MI) {

      const int ind = MI->second;
      M[ind] = MI->first.long0;
      N[ind] = MI->first.long1;
      K[ind] = MI->first.long2;
    }

    /* In order to create a nicer output the profiling data is sorted in
       terms of CPU time spent. Create a vector of pairs for carrying
       out this sort. */
    vector<pair<double, unsigned int> > sortedTime;

    for(unsigned int i=0; i<GEMM_Profile_TotTime.size(); ++i)
      sortedTime.push_back(make_pair(GEMM_Profile_TotTime[i], i));

    sort(sortedTime.begin(), sortedTime.end());

    /* Open the profiling file. */
    ofstream Profile_File;
    Profile_File.precision(15);
    Profile_File.open("gemm_profiling.csv");

    /* Create the CSV header */
    Profile_File << "\"Total_Time\", \"N_Calls\", \"Avg_Time\", \"Min_Time\", \"Max_Time\", \"M\", \"N\", \"K\", \"Avg GFLOPs\"" << endl;

    /* Loop through the different items, where the item with the largest total time is
       written first. As sortedTime is sorted in increasing order, the sequence of
       sortedTime must be reversed. */
    for(vector<pair<double, unsigned int> >::reverse_iterator rit =sortedTime.rbegin();
                                                              rit!=sortedTime.rend(); ++rit) {
      /* Determine the original index in the profiling vectors. */
      const unsigned int ind = rit->second;
      const double AvgTime = GEMM_Profile_TotTime[ind]/GEMM_Profile_NCalls[ind];
      const double GFlops   = 2.0e-9*M[ind]*N[ind]*K[ind]/AvgTime;

      /* Write the data. */
      Profile_File << scientific << GEMM_Profile_TotTime[ind] << ", " << GEMM_Profile_NCalls[ind] << ", "
                   << AvgTime << ", " << GEMM_Profile_MinTime[ind] << ", " << GEMM_Profile_MaxTime[ind] << ", "
                   << M[ind] << ", " << N[ind] << ", " << K[ind] << ", " << GFlops << endl;
    }

    /* Close the file. */
    Profile_File.close();
  }

#endif

}

void CConfig::SetFreeStreamTurboNormal(const su2double* turboNormal){

  FreeStreamTurboNormal[0] = turboNormal[0];
  FreeStreamTurboNormal[1] = turboNormal[1];
  FreeStreamTurboNormal[2] = 0.0;

}

void CConfig::SetMultizone(const CConfig *driver_config, const CConfig* const* config_container){

  for (unsigned short iZone = 0; iZone < nZone; iZone++){

    if (config_container[iZone]->GetTime_Domain() != GetTime_Domain()){
      SU2_MPI::Error("Option TIME_DOMAIN must be the same in all zones.", CURRENT_FUNCTION);
    }
    if (config_container[iZone]->GetnTime_Iter() != GetnTime_Iter()){
      SU2_MPI::Error("Option TIME_ITER must be the same in all zones.", CURRENT_FUNCTION);
    }
    if (config_container[iZone]->GetnOuter_Iter() != GetnOuter_Iter()){
      SU2_MPI::Error("Option OUTER_ITER must be the same in all zones.", CURRENT_FUNCTION);
    }
    if (config_container[iZone]->GetTime_Step() != GetTime_Step()){
      SU2_MPI::Error("Option TIME_STEP must be the same in all zones.", CURRENT_FUNCTION);
    }
    if (config_container[iZone]->GetUnst_CFL() != 0.0){
      SU2_MPI::Error("Option UNST_CFL_NUMBER cannot be used in multizone problems (must be 0),"
                     " use a fixed TIME_STEP instead.", CURRENT_FUNCTION);
    }
    if (config_container[iZone]->GetMultizone_Problem() != GetMultizone_Problem()){
      SU2_MPI::Error("Option MULTIZONE must be the same in all zones.", CURRENT_FUNCTION);
    }
    if (config_container[iZone]->GetMultizone_Mesh() != GetMultizone_Mesh()){
      SU2_MPI::Error("Option MULTIZONE_MESH must be the same in all zones.", CURRENT_FUNCTION);
    }
    if(config_container[iZone]->GetWnd_Cauchy_Crit()){
      SU2_MPI::Error("Option WINDOW_CAUCHY_CRIT must be deactivated for multizone problems.", CURRENT_FUNCTION);
    }
  }
  if(driver_config->GetWnd_Cauchy_Crit()){
    SU2_MPI::Error("Option WINDOW_CAUCHY_CRIT must be deactivated for multizone problems.", CURRENT_FUNCTION);
  }

  bool multiblockDriver = false;
  for (unsigned short iFiles = 0; iFiles < driver_config->GetnVolumeOutputFiles(); iFiles++){
    if (driver_config->GetVolumeOutputFiles()[iFiles] == OUTPUT_TYPE::PARAVIEW_MULTIBLOCK){
      multiblockDriver = true;
    }
  }

  bool multiblockZone = false;
  for (unsigned short iZone = 0; iZone < nZone; iZone++){
    multiblockZone = false;
    for (unsigned short iFiles = 0; iFiles < config_container[iZone]->GetnVolumeOutputFiles(); iFiles++){
      if (config_container[iZone]->GetVolumeOutputFiles()[iFiles] == OUTPUT_TYPE::PARAVIEW_MULTIBLOCK){
        multiblockZone = true;
      }
    }
    if (multiblockZone != multiblockDriver){
      SU2_MPI::Error("To enable PARAVIEW_MULTIBLOCK output, add it to OUTPUT_FILES option in main config and\n"
                     "remove option from sub-config files.", CURRENT_FUNCTION);
    }
  }

  /*--- Fix the Time Step for all subdomains, for the case of time-dependent problems ---*/
  if (driver_config->GetTime_Domain()){
    Delta_UnstTime = driver_config->GetTime_Step();

    Time_Domain = true;
  }

  /*------------------------------------------------------------*/
  /*------ Determine the special properties of the problem -----*/
  /*------------------------------------------------------------*/

  bool fluid_zone = false;
  bool structural_zone = false;

  /*--- If there is at least a fluid and a structural zone ---*/
  for (auto iZone = 0u; iZone < nZone; iZone++) {
    fluid_zone |= config_container[iZone]->GetFluidProblem();
    structural_zone |= config_container[iZone]->GetStructuralProblem();
  }

  if (structural_zone) Relaxation = true;

  /*--- If the problem has FSI properties ---*/
  FSI_Problem = fluid_zone && structural_zone;

  Multizone_Residual = true;
}
