/*!
 * \file ad_structure.hpp
 * \brief Main routines for the algorithmic differentiation (AD) structure.
 * \author T. Albring, J. Blühdorn
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

#include "../code_config.hpp"
#include "../parallelization/omp_structure.hpp"

/*!
 * \namespace AD
 * \brief Contains routines for the reverse mode of AD.
 * In case there is no reverse type configured, they have no effect at all,
 * and so the real versions of the routined are after #else.
 */
namespace AD {
#ifndef CODI_REVERSE_TYPE
/*!
 * \brief Start the recording of the operations and involved variables.
 * If called, the computational graph of all operations occuring after the call will be stored,
 * starting with the variables registered with RegisterInput.
 */
inline void StartRecording() {}

/*!
 * \brief Stops the recording of the operations and variables.
 */
inline void StopRecording() {}

/*!
 * \brief Check if the tape is active
 * \param[out] Boolean which determines whether the tape is active.
 */
inline bool TapeActive() { return false; }

/*!
 * \brief Prints out tape statistics.
 */
inline void PrintStatistics() {}

/*!
 * \brief Registers the variable as an input. I.e. as a leaf of the computational graph.
 * \param[in] data - The variable to be registered as input.
 */
inline void RegisterInput(su2double& data) {}

/*!
 * \brief Registers the variable as an output. I.e. as the root of the computational graph.
 * \param[in] data - The variable to be registered as output.
 */
inline void RegisterOutput(su2double& data) {}

/*!
 * \brief Resize the adjoint vector, for subsequent access without bounds checking.
 */
inline void ResizeAdjoints() {}

/*!
 * \brief Sets the adjoint value at index to val
 * \param[in] index - Position in the adjoint vector.
 * \param[in] val - adjoint value to be set.
 */
inline void SetDerivative(int index, const double val) {}

/*!
 * \brief Extracts the adjoint value at index
 * \param[in] index - position in the adjoint vector where the derivative will be extracted.
 * \return Derivative value.
 */
inline double GetDerivative(int index) { return 0.0; }

/*!
 * \brief Clears the currently stored adjoints but keeps the computational graph.
 */
inline void ClearAdjoints() {}

/*!
 * \brief Computes the adjoints, i.e. the derivatives of the output with respect to the input variables.
 */
inline void ComputeAdjoint() {}

/*!
 * \brief Computes the adjoints, i.e. the derivatives of the output with respect to the input variables.
 * \param[in] enter - Position where we start evaluating the tape.
 * \param[in] leave - Position where we stop evaluating the tape.
 */
inline void ComputeAdjoint(unsigned short enter, unsigned short leave) {}

/*!
 * \brief Computes the adjoints, i.e., the derivatives of the output with respect to the input variables, using forward
 * tape evaluation.
 */
inline void ComputeAdjointForward() {}

/*!
 * \brief Reset the tape structure to be ready for a new recording.
 */
inline void Reset() {}

/*!
 * \brief Reset the variable (set index to zero).
 * \param[in] data - the variable to be unregistered from the tape.
 */
inline void ResetInput(su2double& data) {}

/*!
 * \brief Sets the scalar inputs of a preaccumulation section.
 * \param[in] data - the scalar input variables.
 */
template <class... Ts>
inline void SetPreaccIn(Ts&&... data) {}

/*!
 * \brief Sets the input variables of a preaccumulation section using a 1D array.
 * \param[in] data - the input 1D array.
 * \param[in] size - size of the array.
 */
template <class T>
inline void SetPreaccIn(const T& data, const int size) {}

/*!
 * \brief Sets the input variables of a preaccumulation section using a 2D array.
 * \param[in] data - the input 2D array.
 * \param[in] size_x - size of the array in x dimension.
 * \param[in] size_y - size of the array in y dimension.
 */
template <class T>
inline void SetPreaccIn(const T& data, const int size_x, const int size_y) {}

/*!
 * \brief Starts a new preaccumulation section and sets the input variables.
 *
 * The idea of preaccumulation is to store only the Jacobi matrix of a code section during
 * the taping process instead of all operations. This decreases the tape size and reduces runtime.
 *
 * Input/Output of the section are set with several calls to SetPreaccIn()/SetPreaccOut().
 *
 * Note: the call of this routine must be followed by a call of EndPreacc() and the end of the code section.
 */
inline void StartPreacc() {}

/*!
 * \brief Sets the scalar outputs of a preaccumulation section.
 * \param[in] data - the scalar output variables.
 */
template <class... Ts>
inline void SetPreaccOut(Ts&&... data) {}

/*!
 * \brief Sets the output variables of a preaccumulation section using a 1D array.
 * \param[in] data - the output 1D array.
 */
template <class T>
inline void SetPreaccOut(T&& data, const int size) {}

/*!
 * \brief Sets the input variables of a preaccumulation section using a 2D array.
 * \param[in] data - the output 1D array.
 */
template <class T>
inline void SetPreaccOut(T&& data, const int size_x, const int size_y) {}

/*!
 * \brief Ends a preaccumulation section and computes the local Jacobi matrix
 * of a code section using the variables set with SetLocalInput(), SetLocalOutput() and pushes a statement
 * for each output variable to the AD tape.
 */
inline void EndPreacc() {}

/*!
 * \brief Sets the scalar input of a externally differentiated function.
 * \param[in] data - the scalar input variable.
 */
inline void SetExtFuncIn(su2double& data) {}

/*!
 * \brief Sets the input variables of a externally differentiated function using a 1D array.
 * \param[in] data - the input 1D array.
 * \param[in] size - number of rows.
 */
template <class T>
inline void SetExtFuncIn(const T& data, const int size) {}

/*!
 * \brief  Sets the input variables of a externally differentiated function using a 2D array.
 * \param[in] data - the input 2D array.
 * \param[in] size_x - number of rows.
 * \param[in] size_y - number of columns.
 */
template <class T>
inline void SetExtFuncIn(const T& data, const int size_x, const int size_y) {}

/*!
 * \brief Sets the scalar output of a externally differentiated function.
 * \param[in] data - the scalar output variable.
 */
inline void SetExtFuncOut(su2double& data) {}

/*!
 * \brief Sets the output variables of a externally differentiated function using a 1D array.
 * \param[in] data - the output 1D array.
 * \param[in] size - number of rows.
 */
template <class T>
inline void SetExtFuncOut(T&& data, const int size) {}

/*!
 * \brief  Sets the output variables of a externally differentiated function using a 2D array.
 * \param[in] data - the output 2D array.
 * \param[in] size_x - number of rows.
 * \param[in] size_y - number of columns.
 */
template <class T>
inline void SetExtFuncOut(T&& data, const int size_x, const int size_y) {}

/*!
 * \brief Evaluates and saves gradient data from a variable.
 * \param[in] data - variable whose gradient information will be extracted.
 * \param[in] index - where obtained gradient information will be stored.
 */
inline void SetIndex(int& index, const su2double& data) {}

/*!
 * \brief Pushes back the current tape position to the tape position's vector.
 */
inline void Push_TapePosition() {}

/*!
 * \brief Start a passive region, i.e. stop recording.
 * \return True if tape was active.
 */
inline bool BeginPassive() { return false; }

/*!
 * \brief End a passive region, i.e. start recording if we were recording before.
 * \param[in] wasActive - Whether we were recording before entering the passive region.
 */
inline void EndPassive(bool wasActive) {}

/*!
 * \brief Pause the use of preaccumulation.
 * \return True if preaccumulation was active.
 */
inline bool PausePreaccumulation() { return false; }

/*!
 * \brief Resume the use of preaccumulation.
 * \param[in] wasActive - Whether preaccumulation was active before pausing.
 */
inline void ResumePreaccumulation(bool wasActive) {}

/*!
 * \brief Begin a hybrid parallel adjoint evaluation mode that assumes an inherently safe reverse path.
 */
inline void StartNoSharedReading() {}

/*!
 * \brief End the "no shared reading" adjoint evaluation mode.
 */
inline void EndNoSharedReading() {}

#else
using CheckpointHandler = codi::ExternalFunctionUserData;

using Tape = su2double::Tape;

#ifdef HAVE_OPDI
using ExtFuncHelper = codi::OpenMPExternalFunctionHelper<su2double>;
#else
using ExtFuncHelper = codi::ExternalFunctionHelper<su2double>;
#endif

extern ExtFuncHelper FuncHelper;

extern bool PreaccActive;
#ifdef HAVE_OPDI
SU2_OMP(threadprivate(PreaccActive))
#endif

extern bool PreaccEnabled;

#ifdef HAVE_OPDI
using CoDiTapePosition = Tape::Position;
using OpDiState = void*;
using TapePosition = std::pair<CoDiTapePosition, OpDiState>;
#else
using TapePosition = Tape::Position;
#endif

extern TapePosition StartPosition, EndPosition;

extern std::vector<TapePosition> TapePositions;

extern codi::PreaccumulationHelper<su2double> PreaccHelper;
#ifdef HAVE_OPDI
SU2_OMP(threadprivate(PreaccHelper))
#endif

/*--- Reference to the tape. ---*/

FORCEINLINE Tape& getTape() { return su2double::getTape(); }

FORCEINLINE void RegisterInput(su2double& data) { AD::getTape().registerInput(data); }

FORCEINLINE void RegisterOutput(su2double& data) { AD::getTape().registerOutput(data); }

FORCEINLINE void ResetInput(su2double& data) { data = data.getValue(); }

FORCEINLINE void StartRecording() { AD::getTape().setActive(); }

FORCEINLINE void StopRecording() { AD::getTape().setPassive(); }

FORCEINLINE bool TapeActive() { return AD::getTape().isActive(); }

FORCEINLINE void PrintStatistics() { AD::getTape().printStatistics(); }

FORCEINLINE void ClearAdjoints() { AD::getTape().clearAdjoints(); }

FORCEINLINE void ComputeAdjoint() {
#if defined(HAVE_OPDI)
  opdi::logic->prepareEvaluate();
#endif
  AD::getTape().evaluate();
}

FORCEINLINE void ComputeAdjoint(unsigned short enter, unsigned short leave) {
#if defined(HAVE_OPDI)
  opdi::logic->recoverState(TapePositions[enter].second);
  opdi::logic->prepareEvaluate();
  AD::getTape().evaluate(TapePositions[enter].first, TapePositions[leave].first);
#else
  AD::getTape().evaluate(TapePositions[enter], TapePositions[leave]);
#endif
}

FORCEINLINE void ComputeAdjointForward() { AD::getTape().evaluateForward(); }

FORCEINLINE void Reset() {
  AD::getTape().reset();
#if defined(HAVE_OPDI)
  opdi::logic->reset();
#endif
  if (TapePositions.size() != 0) {
#if defined(HAVE_OPDI)
    for (TapePosition& pos : TapePositions) {
      opdi::logic->freeState(pos.second);
    }
#endif
    TapePositions.clear();
  }
}

FORCEINLINE void ResizeAdjoints() { AD::getTape().resizeAdjointVector(); }

FORCEINLINE void SetIndex(int& index, const su2double& data) { index = data.getIdentifier(); }

// WARNING: For performance reasons, this method does not perform bounds checking.
// When using it, please ensure sufficient adjoint vector size by a call to AD::ResizeAdjoints().
FORCEINLINE void SetDerivative(int index, const double val) {
  if (index == 0)  // Allow multiple threads to "set the derivative" of passive variables without causing data races.
    return;

  AD::getTape().setGradient(index, val, codi::AdjointsManagement::Manual);
}

// WARNING: For performance reasons, this method does not perform bounds checking.
// If called after tape evaluations, the adjoints should exist.
// Otherwise, please ensure sufficient adjoint vector size by a call to AD::ResizeAdjoints().
FORCEINLINE double GetDerivative(int index) {
  return AD::getTape().getGradient(index, codi::AdjointsManagement::Manual);
}

FORCEINLINE bool IsIdentifierActive(su2double const& value) {
  return getTape().isIdentifierActive(value.getIdentifier());
}

/*--- Base case for parameter pack expansion. ---*/
FORCEINLINE void SetPreaccIn() {}

template <class T, class... Ts, su2enable_if<std::is_same<T, su2double>::value> = 0>
FORCEINLINE void SetPreaccIn(const T& data, Ts&&... moreData) {
  if (!PreaccActive) return;
  if (IsIdentifierActive(data)) PreaccHelper.addInput(data);
  SetPreaccIn(moreData...);
}

template <class T, class... Ts, su2enable_if<std::is_same<T, su2double>::value> = 0>
FORCEINLINE void SetPreaccIn(T&& data, Ts&&... moreData) {
  static_assert(!std::is_same<T, su2double>::value, "rvalues cannot be registered");
}

template <class T>
FORCEINLINE void SetPreaccIn(const T& data, const int size) {
  if (PreaccActive) {
    for (int i = 0; i < size; i++) {
      if (IsIdentifierActive(data[i])) {
        PreaccHelper.addInput(data[i]);
      }
    }
  }
}

template <class T>
FORCEINLINE void SetPreaccIn(const T& data, const int size_x, const int size_y) {
  if (!PreaccActive) return;
  for (int i = 0; i < size_x; i++) {
    for (int j = 0; j < size_y; j++) {
      if (IsIdentifierActive(data[i][j])) {
        PreaccHelper.addInput(data[i][j]);
      }
    }
  }
}

FORCEINLINE void StartPreacc() {
  if (AD::getTape().isActive() && PreaccEnabled) {
    PreaccHelper.start();
    PreaccActive = true;
  }
}

/*--- Base case for parameter pack expansion. ---*/
FORCEINLINE void SetPreaccOut() {}

template <class T, class... Ts, su2enable_if<std::is_same<T, su2double>::value> = 0>
FORCEINLINE void SetPreaccOut(T& data, Ts&&... moreData) {
  if (!PreaccActive) return;
  if (IsIdentifierActive(data)) PreaccHelper.addOutput(data);
  SetPreaccOut(moreData...);
}

template <class T>
FORCEINLINE void SetPreaccOut(T&& data, const int size) {
  if (PreaccActive) {
    for (int i = 0; i < size; i++) {
      if (IsIdentifierActive(data[i])) {
        PreaccHelper.addOutput(data[i]);
      }
    }
  }
}

template <class T>
FORCEINLINE void SetPreaccOut(T&& data, const int size_x, const int size_y) {
  if (!PreaccActive) return;
  for (int i = 0; i < size_x; i++) {
    for (int j = 0; j < size_y; j++) {
      if (IsIdentifierActive(data[i][j])) {
        PreaccHelper.addOutput(data[i][j]);
      }
    }
  }
}

FORCEINLINE void Push_TapePosition() {
#if defined(HAVE_OPDI)
  TapePositions.push_back({AD::getTape().getPosition(), opdi::logic->exportState()});
#else
  TapePositions.push_back(AD::getTape().getPosition());
#endif
}

FORCEINLINE void EndPreacc() {
  if (PreaccActive) {
    PreaccHelper.finish(false);
    PreaccActive = false;
  }
}

FORCEINLINE void SetExtFuncIn(const su2double& data) { FuncHelper.addInput(data); }

template <class T>
FORCEINLINE void SetExtFuncIn(const T& data, const int size) {
  for (int i = 0; i < size; i++) {
    FuncHelper.addInput(data[i]);
  }
}

template <class T>
FORCEINLINE void SetExtFuncIn(const T& data, const int size_x, const int size_y) {
  for (int i = 0; i < size_x; i++) {
    for (int j = 0; j < size_y; j++) {
      FuncHelper.addInput(data[i][j]);
    }
  }
}

FORCEINLINE void SetExtFuncOut(su2double& data) {
  if (AD::getTape().isActive()) {
    FuncHelper.addOutput(data);
  }
}

template <class T>
FORCEINLINE void SetExtFuncOut(T&& data, const int size) {
  for (int i = 0; i < size; i++) {
    if (AD::getTape().isActive()) {
      FuncHelper.addOutput(data[i]);
    }
  }
}

template <class T>
FORCEINLINE void SetExtFuncOut(T&& data, const int size_x, const int size_y) {
  for (int i = 0; i < size_x; i++) {
    for (int j = 0; j < size_y; j++) {
      if (AD::getTape().isActive()) {
        FuncHelper.addOutput(data[i][j]);
      }
    }
  }
}

FORCEINLINE void delete_handler(void* handler) {
  CheckpointHandler* checkpoint = static_cast<CheckpointHandler*>(handler);
  checkpoint->clear();
}

FORCEINLINE bool BeginPassive() {
  if (AD::getTape().isActive()) {
    AD::getTape().setPassive();
    return true;
  }
  return false;
}

FORCEINLINE void EndPassive(bool wasActive) {
  if (wasActive) AD::getTape().setActive();
}

FORCEINLINE void StartNoSharedReading() {
#ifdef HAVE_OPDI
  opdi::logic->setAdjointAccessMode(opdi::LogicInterface::AdjointAccessMode::Classical);
  opdi::logic->addReverseBarrier();
#endif
}

FORCEINLINE void EndNoSharedReading() {
#ifdef HAVE_OPDI
  opdi::logic->setAdjointAccessMode(opdi::LogicInterface::AdjointAccessMode::Atomic);
  opdi::logic->addReverseBarrier();
#endif
}

FORCEINLINE bool PausePreaccumulation() {
  const auto current = PreaccEnabled;
  if (!current) return false;
  SU2_OMP_SAFE_GLOBAL_ACCESS(PreaccEnabled = false;)
  return true;
}

FORCEINLINE void ResumePreaccumulation(bool wasActive) {
  if (!wasActive) return;
  SU2_OMP_SAFE_GLOBAL_ACCESS(PreaccEnabled = true;)
}

#endif  // CODI_REVERSE_TYPE

void Initialize();
void Finalize();

}  // namespace AD

/*--- If we compile under OSX we have to overload some of the operators for
 *   complex numbers to avoid the use of the standard operators
 *  (they use a lot of functions that are only defined for doubles) ---*/

#ifdef __APPLE__

namespace std {

template <>
inline su2double abs(const complex<su2double>& x) {
  return sqrt(x.real() * x.real() + x.imag() * x.imag());
}

template <>
inline complex<su2double> operator/(const complex<su2double>& x, const complex<su2double>& y) {
  su2double d = (y.real() * y.real() + y.imag() * y.imag());
  su2double real = (x.real() * y.real() + x.imag() * y.imag()) / d;
  su2double imag = (x.imag() * y.real() - x.real() * y.imag()) / d;

  return complex<su2double>(real, imag);
}

template <>
inline complex<su2double> operator*(const complex<su2double>& x, const complex<su2double>& y) {
  su2double real = (x.real() * y.real() - x.imag() * y.imag());
  su2double imag = (x.imag() * y.real() + x.real() * y.imag());

  return complex<su2double>(real, imag);
}
}  // namespace std
#endif
