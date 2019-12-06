/*!
 * \class CDiscAdjMultizoneDriver.hpp
 * \brief Class for driving adjoint multi-zone problems.
 * \author O. Burghardt, T. Albring, R. Sanchez
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
#include "CMultizoneDriver.hpp"

class CDiscAdjMultizoneDriver : public CMultizoneDriver {

protected:

  /*!
   * \brief Kinds of recordings (three different ones).
   */
  enum class Kind_Tape {
    FULL_TAPE,                /*!< \brief Entire derivative information for a coupled adjoint
                                          solution update. */
    OBJECTIVE_FUNCTION_TAPE,  /*!< \brief Record only the dependence of the objective function
                                          w.r.t. solver variables (from all zones). */
    ZONE_SPECIFIC_TAPE        /*!< \brief Record only the dependence of the solution update in a
                                          specified zone w.r.t. solver variables (from all zones). */
  };

  /*!
   * \brief Position markers within a tape.
   */
  enum Tape_Positions {
    START = 0,                    /*!< \brief Beginning of the tape. */
    REGISTERED = 1,               /*!< \brief Solver variables are registered on the tape. */
    DEPENDENCIES = 2,             /*!< \brief Derived values (e.g. gradients) are set. */
    OBJECTIVE_FUNCTION = 3,       /*!< \brief Objective function is set. */
    TRANSFER = 4,                 /*!< \brief Solution data is transferred between coupled solvers of
                                              different physical zones (giving cross contributions). */
    ITERATION_READY = 5           /*!< \brief Until here, all derivative information is gathered so
                                              that it can be connected to a solver update evaluation. */
  };

  int RecordingState;             /*!< \brief The kind of recording that the tape currently holds. */
  bool retape;                    /*!< \brief Boolean whether a full tape can be kept in memory. */

  su2double ObjFunc;              /*!< \brief Value of the objective function. */
  int ObjFunc_Index;              /*!< \brief Index of the value of the objective function. */

  CIteration*** direct_iteration;       /*!< \brief Array of pointers to the direct iterations. */
  COutput** direct_output;              /*!< \brief Array of pointers to the direct outputs. */
  vector<unsigned short> direct_nInst;  /*!< \brief Total number of instances in the direct problem. */
  vector<unsigned long> nInnerIter;     /*!< \brief Number of inner iterations for each zone. */

  su2vector<bool> Has_Deformation;  /*!< \brief True if iZone has mesh deformation (used for
                                                lazy evaluation of TRANSFER tape section). */

  /*!< \brief Individual cross-terms of the coupled problem, 5D array [iZone][jZone][iSol](iPoint,iVar).
              The column sum, i.e. along all iZones for each jZone, gives the External (total cross-term)
              for jZone, we need to store all terms to have BGS-type updates with relaxation. */
  vector<vector<vector<su2passivematrix> > > Cross_Terms;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] confFile - Configuration file name.
   * \param[in] val_nZone - Total number of zones.
   * \param[in] MPICommunicator - MPI communicator for SU2.
   */
  CDiscAdjMultizoneDriver(char* confFile,
             unsigned short val_nZone,
             SU2_Comm MPICommunicator);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjMultizoneDriver(void);

  /*!
   * \brief [Overload] Launch the computation for discrete adjoint multizone problems.
   */
  void StartSolver() override;

protected:

  /*!
   * \brief [Overload] Run an discrete adjoint update of all solvers within multiple zones.
   */
  void Run() override;

  /*!
   * \brief Evaluate sensitivites for the current adjoint solution and output files.
   * \param[in] iOuterIter - Current outer iteration.
   * \param[in] StopCalc - Final iteration flag (converged or reached max number of iters).
   */
  void EvaluateSensitivities(unsigned long iOuterIter, bool StopCalc);

  /*!
   * \brief Setup the matrix of cross-terms.
   */
  void InitializeCrossTerms();

  /*!
   * \brief Record one iteration of the primal problem within each zone.
   * \param[in] kind_recording - Kind of variables with respect to which we are recording.
   * \param[in] tape_type - indicator which part of a solution update will be recorded.
   * \param[in] record_zone - zone where solution update will be recorded.
   */
  void SetRecording(unsigned short kind_recording, Kind_Tape tape_type, unsigned short record_zone);

  /*!
   * \brief Transfer data between zones and update grids when required.
   */
  void HandleDataTransfer();

  /*!
   * \brief Run one direct iteration in a zone.
   * \param[in] iZone - Zone in which we run an iteration.
   * \param[in] kind_recording - Kind of variables with respect to which we are recording.
   */
  void DirectIteration(unsigned short iZone, unsigned short kind_recording);

  /*!
   * \brief Set the objective function.
   * \param[in] kind_recording - Kind of variables with respect to which we are recording.
   */
  void SetObjFunction(unsigned short kind_recording);

  /*!
   * \brief Initialize the adjoint value of the objective function.
   */
  void SetAdj_ObjFunction();

  /*!
   * \brief Summary of all routines to evaluate the adjoints in iZone.
   * \param[in] iZone - Zone in which adjoints are evaluated depending on their (preceding) seeding.
   * \param[in] eval_transfer - Evaluate adjoints of transfer and mesh deformation routines.
   */
  void ComputeAdjoints(unsigned short iZone, bool eval_transfer = true);

  /*!
   * \brief Puts BGSSolution_k back into Solution.
   * \param[in] iZone - Zone index.
   */
  void Set_Solution_To_BGSSolution_k(unsigned short iZone);

  /*!
   * \brief Puts Solution into BGSSolution_k.
   * \param[in] iZone - Zone index.
   */
  void Set_BGSSolution_k_To_Solution(unsigned short iZone);

  /*!
   * \brief Add Solution vector to External.
   * \param[in] iZone - Zone index.
   */
  void Add_Solution_To_External(unsigned short iZone);

  /*!
   * \brief Add External_Old vector to Solution.
   * \param[in] iZone - Zone index.
   */
  void Add_External_To_Solution(unsigned short iZone);

  /*!
   * \brief Extract contribution of iZone to jZone with BGS relaxation.
   * \param[in] iZone - Source zone (the one that was initialized).
   * \param[in] jZone - Target zone (the one that transfers to iZone in the primal problem).
   */
  void Update_Cross_Term(unsigned short iZone, unsigned short jZone);

  /*!
   * \brief Saves the current (adjoint) Solution vector to Solution_BGS_k.
   * \param[in] iZone - Zone index.
   */
  void Set_BGSSolution(unsigned short iZone);

  /*!
   * \brief Compute BGS residuals.
   * \param[in] iZone - Zone where solver residuals are computed.
   */
  void SetResidual_BGS(unsigned short iZone);
};
