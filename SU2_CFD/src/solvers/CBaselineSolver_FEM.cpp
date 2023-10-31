/*!
 * \file CBaselineSolver_FEM.cpp
 * \brief Main subroutines for CBaselineSolver_FEM class.
 * \author F. Palacios, T. Economon
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


#include "../../include/solvers/CBaselineSolver_FEM.hpp"


CBaselineSolver_FEM::CBaselineSolver_FEM() : CSolver() { }

CBaselineSolver_FEM::CBaselineSolver_FEM(CGeometry *geometry, CConfig *config) {

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. If necessary, it is
   possible to increase nMatchingFacesWithHaloElem a bit, such that
   the computation of the external faces may be more efficient when
   using multiple threads. ---*/

  auto *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  nVolElemTot   = DGGeometry->GetNVolElemTot();
  nVolElemOwned = DGGeometry->GetNVolElemOwned();
  volElem       = DGGeometry->GetVolElem();

  /*--- Routines to access the number of variables and string names. ---*/

  SetOutputVariables(geometry, config);

  /*--- Determine the total number of DOFs stored on this rank and allocate the memory
   to store the conservative variables. ---*/
  nDOFsLocOwned = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) nDOFsLocOwned += volElem[i].nDOFsSol;

  nDOFsLocTot = nDOFsLocOwned;
  for(unsigned long i=nVolElemOwned; i<nVolElemTot; ++i) nDOFsLocTot += volElem[i].nDOFsSol;

  VecSolDOFs.resize(nVar*nDOFsLocTot);

  /*--- Determine the global number of DOFs. ---*/
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nDOFsLocOwned, &nDOFsGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
#else
  nDOFsGlobal = nDOFsLocOwned;
#endif

  /*--- Store the number of DOFs in the geometry class in case of restart. ---*/
  geometry->SetnPointDomain(nDOFsLocOwned);
  geometry->SetGlobal_nPointDomain(nDOFsGlobal);

  /*--- Initialize the solution to zero. ---*/

  unsigned long ii = 0;
  for(unsigned long i=0; i<nDOFsLocTot; ++i) {
    for(unsigned short j=0; j<nVar; ++j, ++ii) {
      VecSolDOFs[ii] = 0.0;
    }
  }

}

void CBaselineSolver_FEM::SetOutputVariables(CGeometry *geometry, CConfig *config) {

  /*--- Open the restart file and extract the nVar and field names. ---*/

  string Tag, text_line, AdjExt, UnstExt;
  unsigned long TimeIter = config->GetTimeIter();

  ifstream restart_file;
  string filename;

  /*--- Retrieve filename from config ---*/

  filename = config->GetSolution_FileName();

  /*--- Unsteady problems require an iteration number to be appended. ---*/

  if (config->GetTime_Domain()) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(TimeIter), ".dat");
  }

  /*--- Read only the number of variables in the restart file. ---*/

  if (config->GetRead_Binary_Restart()) {

    int nVar_Buf = 5;
    int var_buf[5];

#ifndef HAVE_MPI

    /*--- Serial binary input. ---*/

    FILE *fhw;
    fhw = fopen(filename.c_str(),"rb");
    size_t ret;

    /*--- Error check for opening the file. ---*/

    if (!fhw)
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + filename,
                     CURRENT_FUNCTION);

    /*--- First, read the number of variables and points. ---*/

    ret = fread(var_buf, sizeof(int), nVar_Buf, fhw);
    if (ret != (unsigned long)nVar_Buf) {
      SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
    }

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (var_buf[0] != 535532)
      SU2_MPI::Error(string("File ") + filename + string(" is not a binary SU2 restart file.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the READ_BINARY_RESTART option."), CURRENT_FUNCTION);

    /*--- Close the file. ---*/

    fclose(fhw);

#else

    /*--- Parallel binary input using MPI I/O. ---*/

    MPI_File fhw;
    int ierr;

    /*--- All ranks open the file using MPI. ---*/

    char fname[100];
    strcpy(fname, filename.c_str());
    ierr = MPI_File_open(SU2_MPI::GetComm(), fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

    /*--- Error check opening the file. ---*/

    if (ierr)
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + filename,
                     CURRENT_FUNCTION);

    /*--- First, read the number of variables and points (i.e., cols and rows),
     which we will need in order to read the file later. Also, read the
     variable string names here. Only the master rank reads the header. ---*/

    if (rank == MASTER_NODE)
      MPI_File_read(fhw, var_buf, nVar_Buf, MPI_INT, MPI_STATUS_IGNORE);

    /*--- Broadcast the number of variables to all procs and store more clearly. ---*/

    SU2_MPI::Bcast(var_buf, nVar_Buf, MPI_INT, MASTER_NODE, SU2_MPI::GetComm());

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (var_buf[0] != 535532)
      SU2_MPI::Error(string("File ") + filename + string(" is not a binary SU2 restart file.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the READ_BINARY_RESTART option."), CURRENT_FUNCTION);

    /*--- All ranks close the file after writing. ---*/

    MPI_File_close(&fhw);

#endif

    /*--- Set the number of variables, one per field in the
     restart file (without including the PointID) ---*/

    nVar = var_buf[1];

  } else {

    /*--- First, check that this is not a binary restart file. ---*/

    int magic_number;

#ifndef HAVE_MPI

    /*--- Serial binary input. ---*/

    FILE *fhw;
    fhw = fopen(filename.c_str(), "rb");
    size_t ret;

    /*--- Error check for opening the file. ---*/

    if (!fhw)
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + filename,
                     CURRENT_FUNCTION);

    /*--- Attempt to read the first int, which should be our magic number. ---*/

    ret = fread(&magic_number, sizeof(int), 1, fhw);
    if (ret != 1) {
      SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
    }

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (magic_number == 535532)
      SU2_MPI::Error(string("File ") + filename + string(" is a binary SU2 restart file, expected ASCII.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the READ_BINARY_RESTART option."), CURRENT_FUNCTION);
    fclose(fhw);

#else

    /*--- Parallel binary input using MPI I/O. ---*/

    MPI_File fhw;
    int ierr;

    /*--- All ranks open the file using MPI. ---*/

    char fname[100];
    strcpy(fname, filename.c_str());
    ierr = MPI_File_open(SU2_MPI::GetComm(), fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

    /*--- Error check opening the file. ---*/

    if (ierr)
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + filename,
                     CURRENT_FUNCTION);

    /*--- Have the master attempt to read the magic number. ---*/

    if (rank == MASTER_NODE)
      MPI_File_read(fhw, &magic_number, 1, MPI_INT, MPI_STATUS_IGNORE);

    /*--- Broadcast the number of variables to all procs and store clearly. ---*/

    SU2_MPI::Bcast(&magic_number, 1, MPI_INT, MASTER_NODE, SU2_MPI::GetComm());

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (magic_number == 535532)
      SU2_MPI::Error(string("File ") + filename + string(" is a binary SU2 restart file, expected ASCII.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the READ_BINARY_RESTART option."), CURRENT_FUNCTION);

    MPI_File_close(&fhw);

#endif

    /*--- Open the restart file ---*/

    restart_file.open(filename.data(), ios::in);

    /*--- In case there is no restart file ---*/

    if (restart_file.fail())
      SU2_MPI::Error(string("SU2 solution file ") + filename + string(" not found"), CURRENT_FUNCTION);

    /*--- Identify the number of fields (and names) in the restart file ---*/

    getline (restart_file, text_line);

    stringstream ss(text_line);
    while (ss >> Tag) {
      config->fields.push_back(Tag);
      if (ss.peek() == ',') ss.ignore();
    }

    /*--- Close the file (the solution date is read later). ---*/

    restart_file.close();

    /*--- Set the number of variables, one per field in the
     restart file (without including the PointID) ---*/

    nVar = config->fields.size() - 1;

    /*--- Clear the fields vector since we'll read it again. ---*/

    config->fields.clear();

  }

}

void CBaselineSolver_FEM::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/
  unsigned short iVar;
  unsigned long index;

  string restart_filename = config->GetSolution_FileName();

  if (config->GetTime_Domain()) {
    restart_filename = config->GetUnsteady_FileName(restart_filename, SU2_TYPE::Int(val_iter), "");
  }

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned short rbuf_NotMatching = 0;
  unsigned long nDOF_Read = 0;

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Load data from the restart into correct containers. ---*/

  counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1];
      for (iVar = 0; iVar < nVar; iVar++) {
        VecSolDOFs[nVar*iPoint_Local+iVar] = Restart_Data[index+iVar];
      }
      /*--- Update the local counter nDOF_Read. ---*/
      ++nDOF_Read;

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- Detect a wrong solution file ---*/
  if(nDOF_Read < nDOFsLocOwned) rbuf_NotMatching = 1;

#ifdef HAVE_MPI
  unsigned short sbuf_NotMatching = rbuf_NotMatching;
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_MAX, SU2_MPI::GetComm());
#endif

  if (rbuf_NotMatching != 0)
    SU2_MPI::Error(string("The solution file ") + restart_filename +
                   string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."),
                   CURRENT_FUNCTION);

  /*--- Delete the class memory that is used to load the restart. ---*/

  delete [] Restart_Vars;
  delete [] Restart_Data;
  Restart_Vars = nullptr; Restart_Data = nullptr;

}

CBaselineSolver_FEM::~CBaselineSolver_FEM() = default;
