/*!
 * \file output_driver.cpp
 * \brief Main subroutines for multizone output
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

#include "../include/output_driver.hpp"

CDriverOutput::CDriverOutput(CConfig **config) {

  nZone = config[ZONE_0]->GetnZone();

  OuterConvergenceTable = new PrintingToolbox::CTablePrinter(&std::cout);

  nRequestedScreenFields = 0;
  nRequestedHistoryFields = 0;

  field_width = 12;

}

CDriverOutput::~CDriverOutput() {


}

void CDriverOutput::SetHeader(COutput **output, CSolver *****solver, CConfig **config) {


}

void CDriverOutput::SetBody(COutput **output, CSolver *****solver, CConfig *driver_config, CConfig **config) {

  /*--- Output using only the master node ---*/

  if (rank == MASTER_NODE) {

    bool write_header, write_history, write_screen;

    /*--- Retrieve residual and extra data -----------------------------------------------------------------*/

    LoadHistoryData(output, solver, driver_config);

    /*--- Write the history file ---------------------------------------------------------------------------*/
    //write_history = WriteHistoryFile_Output(config[val_iZone], DualTime_Iteration);
    //if (write_history) SetHistoryFile_Output(config[val_iZone]);

    /*--- Write the screen header---------------------------------------------------------------------------*/
    write_header = true;
    if (write_header) SetScreen_Header(driver_config, config);

    /*--- Write the screen output---------------------------------------------------------------------------*/
    write_screen = true;
    if (write_screen) SetScreen_Output(output, driver_config, config);

  }


}

void CDriverOutput::LoadHistoryData(COutput **output, CSolver *****solver, CConfig *driver_config) {

  unsigned short iZone, iSol, iVar;
  unsigned short nVar, nVar_Zone;
  string name, header;
  su2double val, avgsol, avgzone;

  SetHistoryOutputValue("TIME_ITER", driver_config->GetTimeIter());
  SetHistoryOutputValue("OUTER_ITER", driver_config->GetOuterIter());

  for (iZone = 0; iZone < nZone; iZone++){

    /*--- Initialize values per zone ---*/
    nVar_Zone = 0;
    avgzone = 0.0;

    /*--- Accounting for all the solvers ---*/
    for (iSol = 0; iSol < MAX_SOLS; iSol++){
      /*-- If the solver position iSol is enabled --*/
      if (solver[iZone][INST_0][MESH_0][iSol] != NULL){

        /*---Initialize values per solver ---*/
        nVar = solver[iZone][INST_0][MESH_0][iSol]->GetnVar();
        avgsol = 0.0;

        /*-- For all the variables per solver --*/
        for (iVar = 0; iVar < nVar; iVar++){

          /*--- Get the unique name for the history data ---*/
          name = "ZONE" + to_string(iZone) + "_SOL" + to_string(iSol) + "_VAL" + to_string(iVar);

          /*--- Set the value of the BGS residual for the variable iVar, solver iSol, zone iZone ---*/
          val = log10(solver[iZone][INST_0][MESH_0][iSol]->GetRes_BGS(iVar));
          SetHistoryOutputValue(name, val);

          /*--- Add values to averages ---*/
          avgsol += val;
          avgzone += val;
          nVar_Zone++;
        }

        /*--- Get the unique name for the averaged history data per solver ---*/
        name = "ZONE" + to_string(iZone) + "_SOL" + to_string(iSol);

        /*--- Compute the average and set the value for the solver iSol, zone iZone---*/
        avgsol = avgsol / nVar;
        SetHistoryOutputValue(name, avgsol);

      }

    }

    /*--- Get an unique name for the averaged history data per zone ---*/
    name = "ZONE" + to_string(iZone);

    /*--- Compute the average and set the value for the zone iZone---*/
    avgzone = avgzone / nVar_Zone;
    SetHistoryOutputValue(name, avgzone);

  }


}

void CDriverOutput::SetHistoryOutputFields(COutput **output, CSolver *****solver, CConfig **config) {

  unsigned short iZone, iSol, iVar;
  unsigned short iReqField;
  string name, header;

  AddHistoryOutput("TIME_ITER", "Time_Iter", FORMAT_INTEGER,  "ITER");
  RequestedScreenFields.push_back("TIME_ITER"); RequestedHistoryFields.push_back("TIME_ITER");
  AddHistoryOutput("OUTER_ITER", "Outer_Iter", FORMAT_INTEGER,  "ITER");
  RequestedScreenFields.push_back("OUTER_ITER"); RequestedHistoryFields.push_back("OUTER_ITER");

  /*--- Set the fields ---*/
  for (iZone = 0; iZone < nZone; iZone++){

    /*--- Accounting for all the solvers ---*/
    for (iSol = 0; iSol < MAX_SOLS; iSol++){
      /*-- If the solver position iSol is enabled --*/
      if (solver[iZone][INST_0][MESH_0][iSol] != NULL){
        /*-- For all the variables per solver --*/
        for (iVar = 0; iVar < solver[iZone][INST_0][MESH_0][iSol]->GetnVar(); iVar++){

          /*--- Set an unique name for the history data ---*/
          name = "ZONE" + to_string(iZone) + "_SOL" + to_string(iSol) + "_VAL" + to_string(iVar);
          /*--- Set an unique name for the history headers ---*/
          header = "res[" + to_string(iZone) + "][" + to_string(iSol) + "][" + to_string(iVar) + "]";

          AddHistoryOutput(name, header, FORMAT_FIXED,  "VAR_RES", TYPE_RESIDUAL);

          /*--- Request the variable residual for history output ---*/
          RequestedHistoryFields.push_back(name);
        }

        /*--- Set an unique name for the averaged history data per solver ---*/
        name = "ZONE" + to_string(iZone) + "_SOL" + to_string(iSol);
        /*--- Set an unique name for the history headers of the averaged data per solver ---*/
        header = solver[iZone][INST_0][MESH_0][iSol]->GetSolverName();
        header += "[" + to_string(iZone) + "]";

        /*--- Add the average residual of the current solver to output ---*/
        AddHistoryOutput(name, header, FORMAT_FIXED,  "SOL_AVGRES", TYPE_RESIDUAL);

        /*--- Request the average residual for screen output ---*/
        RequestedScreenFields.push_back(name);

      }
    }

    /*--- Set an unique name for the averaged history data ---*/
    name = "ZONE" + to_string(iZone);
    /*--- Set an unique name for the history headers of the averaged data ---*/
    header = "avgres[" + to_string(iZone) + "]";

    AddHistoryOutput(name, header, FORMAT_FIXED,  "ZONE_AVGRES", TYPE_RESIDUAL);
  }

  /*--- Retrieve number of requested fields ---*/
  nRequestedScreenFields = RequestedScreenFields.size();
  nRequestedHistoryFields = RequestedHistoryFields.size();

  string RequestedField;

  /*--- Set screen convergence output header ---*/
  for (iReqField = 0; iReqField < nRequestedScreenFields; iReqField++){
    RequestedField = RequestedScreenFields[iReqField];
    if (HistoryOutput_Map.count(RequestedField) > 0){
      OuterConvergenceTable->AddColumn(HistoryOutput_Map[RequestedField].FieldName, field_width);
    } else {
      SU2_MPI::Error(string("Requested screen output field not found: ") + RequestedField, CURRENT_FUNCTION);
    }
  }


}

void CDriverOutput::SetScreen_Header(CConfig *driver_config, CConfig **config) {

  OuterConvergenceTable->PrintHeader();

}


void CDriverOutput::SetScreen_Output(COutput **output, CConfig *driver_config, CConfig **config) {

  string RequestedField;

  for (unsigned short iReqField = 0; iReqField < nRequestedScreenFields; iReqField++){
    stringstream out;
    RequestedField = RequestedScreenFields[iReqField];
    if (HistoryOutput_Map.count(RequestedField) > 0){
      switch (HistoryOutput_Map[RequestedField].ScreenFormat) {
        case FORMAT_INTEGER:
          output[ZONE_0]->PrintScreenInteger(out, SU2_TYPE::Int(HistoryOutput_Map[RequestedField].Value));
          break;
        case FORMAT_FIXED:
          output[ZONE_0]->PrintScreenFixed(out, HistoryOutput_Map[RequestedField].Value);
          break;
        case FORMAT_SCIENTIFIC:
          output[ZONE_0]->PrintScreenScientific(out, HistoryOutput_Map[RequestedField].Value);
          break;
      }
    }
    (*OuterConvergenceTable) << out.str();
  }
}
