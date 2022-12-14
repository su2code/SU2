/*!
 * \file CDiscAdjSinglezoneDriver.hpp
 * \brief Headers of the main subroutines for driving single or multi-zone problems.
 *        The subroutines and functions are in the <i>driver_structure.cpp</i> file.
 * \author T. Economon, H. Kline, R. Sanchez
 * \version 7.4.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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
#include "CSinglezoneDriver.hpp"
#include "../../../Common/include/toolboxes/CQuasiNewtonInvLeastSquares.hpp"
#include "../../../Common/include/linear_algebra/CPreconditioner.hpp"
#include "../../../Common/include/linear_algebra/CMatrixVectorProduct.hpp"
#include "../../../Common/include/linear_algebra/CSysSolve.hpp"

/*!
 * \class CDiscAdjSinglezoneDriver
 * \ingroup DiscAdj
 * \brief Class for driving single-zone adjoint solvers.
 * \author R. Sanchez
 * \version 7.4.0 "Blackbird"
 */
class CDiscAdjSinglezoneDriver : public CSinglezoneDriver {
protected:
#ifdef CODI_FORWARD_TYPE
    using Scalar = su2double;
#else
    using Scalar = passivedouble;
#endif
    
    unsigned long nAdjoint_Iter;                  /*!< \brief The number of adjoint iterations that are run on the fixed-point solver.*/
    RECORDING RecordingState;                     /*!< \brief The kind of recording the tape currently holds.*/
    RECORDING MainVariables;                      /*!< \brief The kind of recording linked to the main variables of the problem.*/
    RECORDING SecondaryVariables;                 /*!< \brief The kind of recording linked to the secondary variables of the problem.*/
    int MainSolver;                               /*!< \brief Index of the main adjoint solver. */
    su2double ObjFunc;                            /*!< \brief The value of the objective function.*/
    CIteration* direct_iteration;                 /*!< \brief A pointer to the direct iteration.*/
    
    CConfig *config;                              /*!< \brief Definition of the particular problem. */
    CIteration *iteration;                        /*!< \brief Container vector with all the iteration methods. */
    CIntegration **integration;                   /*!< \brief Container vector with all the integration methods. */
    CGeometry *geometry;                          /*!< \brief Geometrical definition of the problem. */
    CSolver **solver;                             /*!< \brief Container vector with all the solutions. */
    COutput *direct_output;
    CNumerics ***numerics;                        /*!< \brief Container vector with all the numerics. */
    
    COutputLegacy* output_legacy;
    
    /*!< \brief Fixed-Point corrector that can be applied to inner iterations of the residual-based formulation. */
    CQuasiNewtonInvLeastSquares<passivedouble> Corrector;

    /*!< \brief Members to use GMRES to drive inner iterations (alternative to quasi-Newton) of the residual-based formulation. */
    static constexpr unsigned long KrylovMinIters = 3;
    const Scalar KrylovSysTol = 0.00001;
    const Scalar KrylovPreTol = 0.1;
    bool KrylovSet = false;

    CSysMatrix<Scalar> CopiedJacobian; 
    CSysSolve<Scalar> AdjSolver;
    CSysVector<Scalar> AdjRHS;
    CSysVector<Scalar> AdjSol;
    CPreconditioner<Scalar>* PrimalPreconditioner = nullptr;
    CSysMatrixVectorProduct<Scalar>* PrimalJacobian = nullptr;

    class LinOperator : public CMatrixVectorProduct<Scalar> {
    public:
        CDiscAdjSinglezoneDriver* const driver;

        LinOperator(CDiscAdjSinglezoneDriver* d) : driver(d) { }

        inline void operator()(const CSysVector<Scalar> & u, CSysVector<Scalar> & v) const override {
            driver->ApplyOperator(u, v);
        }
    };

    class LinPreconditioner : public CPreconditioner<Scalar> {
    public:    
        CDiscAdjSinglezoneDriver* const driver;
        
        LinPreconditioner(CDiscAdjSinglezoneDriver* d): driver(d) { }
            
        inline void operator()(const CSysVector<Scalar> & u, CSysVector<Scalar> & v) const override {
            driver->ApplyPreconditioner(u, v);
        }
    };
    
public:
    
    /*!
     * \brief Constructor of the class.
     * \param[in] confFile - Configuration file name.
     * \param[in] val_nZone - Total number of zones.
     * \param[in] val_nDim - Total number of dimensions.
     * \param[in] MPICommunicator - MPI communicator for SU2.
     */
    CDiscAdjSinglezoneDriver(char* confFile,
                             unsigned short val_nZone,
                             SU2_Comm MPICommunicator);
    
    /*!
     * \brief Destructor of the class.
     */
    ~CDiscAdjSinglezoneDriver(void) override;
    
    /*!
     * \brief Preprocess the single-zone iteration
     * \param[in] TimeIter - index of the current time-step.
     */
    void Preprocess(unsigned long TimeIter) override;
    
    /*!
     * \brief Run the discrete adjoint solver with a single zone.
     */
    void Run(void) override;
    
    /*!
     * \brief Update the discrete adjoint solution.
     */
    void UpdateAdjoints(void);

    /*!
     * \brief Update the primal time iteration.
     */
    void UpdateTimeIter(void);
    
    /*!
     * \brief Update the primal far-field variables.
     */
    void UpdateFarfield(void);
    
    /*!
     * \brief Update the primal geometry (does not include mesh deformation).
     */
    void UpdateGeometry(void);

    /*!
     * \brief Deform the primal mesh.
     */
    void DeformGeometry(void);
    
    /*!
     * \brief Update the primal states.
     */
    void UpdateStates();
    
    /*!
     * \brief Update the primal residuals.
     */
    void UpdateResiduals();
    
    /*!
     * \brief Update the primal tractions.
     */
    void UpdateTractions();

    /*!
     * \brief Update the primal objective.
     */
    void UpdateObjective();

    /*!
     * \brief Update the primal jacobian.
     */
    void UpdateJacobians();
    
    /*!
     * \brief Postprocess the adjoint iteration for ZONE_0.
     */
    void Postprocess(void) override;
    
    /*!
     * \brief Record one iteration of a flow iteration in within multiple zones.
     * \param[in] kind_recording - Type of recording (full list in ENUM_RECORDING, option_structure.hpp)
     */
    void SetRecording(RECORDING kind_recording);
    
    /*!
     * \brief Run one iteration of the solver.
     * \param[in] kind_recording - Type of recording (full list in ENUM_RECORDING, option_structure.hpp)
     */
    void DirectRun_FixedPoint(RECORDING kind_recording);
    
    /*!
     * \brief Run one iteration of the solver.
     * \param[in] kind_recording - Type of recording (full list in ENUM_RECORDING, option_structure.hpp)
     */
    void DirectRun_Residual(RECORDING kind_recording);
    
    /*!
     * \brief Print out the direct residuals.
     * \param[in] kind_recording - Type of recording (full list in ENUM_RECORDING, option_structure.hpp)
     */
    void Print_DirectResidual(RECORDING kind_recording);
    
    /*!
     * \brief Record the main computational path.
     */
    void MainRecording(void);
    
    /*!
     * \brief Record the secondary computational path.
     */
    void SecondaryRecording(void);
    
    /*!
     * \brief gets Convergence on physical time scale, (deactivated in adjoint case)
     * \return false
     */
    inline bool GetTimeConvergence() const override {return false;};
    
    /*!
     * \brief Sum the number adjoint variables for all solvers.
     * \return Total number of solution variables.
     */
    unsigned short GetTotalNumberOfVariables() const {
        unsigned short nVar = 0;

        for (auto iSol = 0u; iSol < MAX_SOLS; iSol++) {
            auto solver = solver_container[ZONE_0][INST_0][MESH_0][iSol];

            if (!solver || !solver->GetAdjoint()) continue;
            
            nVar += solver->GetnVar();
        }

        return nVar;
    }

    /*!
     * \brief Get the adjoint states from all solvers.
     * \param[out] values - Values object with interface (iPoint, iVar).
     */
    template<class Container>
    void GetAllSolutions(Container& values) const {
        const auto nPoint = geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();

        /*--- Get all the adjoint state variables ---*/
        unsigned short offset = 0;

        for (auto iSol = 0u; iSol < MAX_SOLS; ++iSol) {
            auto solver = solver_container[ZONE_0][INST_0][MESH_0][iSol];
        
            if (!solver || !solver->GetAdjoint()) continue;
            
            for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
                for (auto iVar = 0ul; iVar < solver->GetnVar(); ++iVar) {
                    auto value = solver->GetNodes()->GetSolution(iPoint, iVar);

                    values(iPoint, offset + iVar) = SU2_TYPE::GetValue(value);
                }
            }

            offset += solver->GetnVar();
        }
    }

    /*!
     * \brief Set the adjoint states from all solvers.
     * \param[out] values - Values object with interface (iPoint, iVar).
     */
    template<class Container>
    void SetAllSolutions(Container& values) const {
        const auto nPoint = geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();

        /*--- Set all the adjoint state variables ---*/
        unsigned short offset = 0;

        for (auto iSol = 0u; iSol < MAX_SOLS; ++iSol) {
            auto solver = solver_container[ZONE_0][INST_0][MESH_0][iSol];
        
            if (!solver || !solver->GetAdjoint()) continue;
            
            for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
                for (auto iVar = 0ul; iVar < solver->GetnVar(); ++iVar) {
                    auto value = values(iPoint, offset + iVar);
                    
                    solver->GetNodes()->SetSolution(iPoint, iVar, value);
                }
            }

            offset += solver->GetnVar();
        }
    }

    /*!
     * \brief Get the partial objective sensitivities from all solvers.
     * \param[out] values - Values object with interface (iPoint, iVar).
     */
    template<class Container>
    void GetAllObjectiveStatesSensitivities(Container& values) const {
        const auto nPoint = geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();

        /*--- Get all the partial sensitivities (dObjective/dStates) ---*/
        unsigned short offset = 0;

        for (auto iSol = 0u; iSol < MAX_SOLS; ++iSol) {
            auto solver = solver_container[ZONE_0][INST_0][MESH_0][iSol];
        
            if (!solver || !solver->GetAdjoint()) continue;
            
            for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
                for (auto iVar = 0ul; iVar < solver->GetnVar(); ++iVar) {
                    auto value = solver->GetSens_dObjective_dStates(iPoint, iVar);

                    values(iPoint, offset + iVar) = -SU2_TYPE::GetValue(value);
                }
            }

            offset += solver->GetnVar();
        }
    }

    /*!
     * \brief Get the partial jacobian-adjoint products from all solvers.
     * \param[out] values - Values object with interface (iPoint, iVar).
     */
    template<class Container>
    void GetAllResidualsStatesSensitivities(Container& values) const {
        const auto nPoint = geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint();

        /*--- Get all the partial jacobian-adjoint products (dResiduals/dStates) ---*/
        unsigned short offset = 0;

        for (auto iSol = 0u; iSol < MAX_SOLS; ++iSol) {
            auto solver = solver_container[ZONE_0][INST_0][MESH_0][iSol];
        
            if (!solver || !solver->GetAdjoint()) continue;
            
            for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
                for (auto iVar = 0ul; iVar < solver->GetnVar(); ++iVar) {
                    auto value = solver->GetProd_dResiduals_dStates(iPoint, iVar);

                    values(iPoint, offset + iVar) = SU2_TYPE::GetValue(value);
                }
            }

            offset += solver->GetnVar();
        }
    }
    
protected:
    
    /*!
     * \brief Run a single iteration of the main fixed-point discrete adjoint solver with a single zone.
     */
    void Run_FixedPoint(void);
    
    /*!
     * \brief Run the computation of the main residual-based discrete adjoint sensitivities with a single zone.
     */
    void Run_Residual(void);
    
    /*!
     * \brief Run a single iteration of the secondary fixed-point discrete adjoint solver with a single zone.
     */
    void SecondaryRun_FixedPoint(void);
    
    /*!
     * \brief Run the computation of the secondary residual-based discrete adjoint sensitivities with a single zone.
     */
    void SecondaryRun_Residual(void);
    
    /*!
     * \brief Update the fixed-point discrete adjoint solver with a single zone.
     */
    void UpdateAdjoints_FixedPoint(void);
    
    /*!
     * \brief Update the residual-based discrete adjoint solver with a single zone.
     */
    void UpdateAdjoints_Residual(void);
    
    /*!
     * \brief Update the residual-based discrete adjoint solver with a single zone.
     */
    void UpdateAdjoints_Objective(void);

    /*!
    * \brief Adjoint problem Jacobian-vector product.
    */
    void ApplyOperator(const CSysVector<Scalar>& u, CSysVector<Scalar>& v);

    /*!
    * \brief Adjoint problem preconditioner (based on the transpose approximate Jacobian of the primal problem).
    */
    void ApplyPreconditioner(const CSysVector<Scalar>& u, CSysVector<Scalar>& v);
};
