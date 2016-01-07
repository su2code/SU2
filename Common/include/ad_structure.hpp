#pragma once

#include "../include/datatype_structure.hpp"

/*!
 * \namespace AD
 * \brief Contains routines for the reverse mode of AD.
 * In case there is no reverse type configured, they have no effect at all.
 */

namespace AD{
  /*!
   * \brief Start the recording of the operations and involved variables.
   * If called, the computational graph of all operations occuring after the call will be stored,
   * starting with the variables registered with RegisterInput.
   */
  void StartRecording();

  /*!
   * \brief Stops the recording of the operations and variables.
   */
  void StopRecording();

  /*!
   * \brief Registers the variable as an input. I.e. as a leaf of the computational graph.
   * \param[in] data - The variable to be registered as input.
   */
  void RegisterInput(su2double &data);

  /*!
   * \brief Registers the variable as an output. I.e. as the root of the computational graph.
   * \param[in] data - The variable to be registered as output.
   */
  void RegisterOutput(su2double &data);

  /*!
   * \brief Clears the currently stored adjoints but keeps the computational graph.
   */
  void ClearAdjoints();

  /*!
   * \brief Computes the adjoints, i.e. the derivatives of the output with respect to the input variables.
   */
  void ComputeAdjoint();

  /*!
   * \brief Reset the tape structure to be ready for a new recording.
   */
  void Reset();

  /*!
   * \brief Reset the variable (set index to zero).
   * \param[in] data - the variable to be unregistered from the tape.
   */
  void ResetInput(su2double &data);

  /*!
   * \brief Sets the scalar input of a preaccumulation section.
   * \param[in] data - the scalar input variable.
   */
  void SetPreaccIn(const su2double &data);

  /*!
   * \brief Sets the input variables of a preaccumulation section using a 1D array.
   * \param[in] data - the input 1D array.
   * \param[in] size - size of the array.
   */
  void SetPreaccIn(su2double* data, const int size);

  /*!
   * \brief Sets the input variables of a preaccumulation section using a 2D array.
   * \param[in] data - the input 2D array.
   * \param[in] size_x - size of the array in x dimension.
   * \param[in] size_y - size of the array in y dimension.
   */
  void SetPreaccIn(su2double** data, const int size_x, const int size_y);

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
  void StartPreacc();

  /*!
   * \brief Sets the scalar output of a preaccumulation section.
   * \param[in] data - the scalar output variable.
   */
  void SetPreaccOut(su2double &data);

  /*!
   * \brief Sets the output variables of a preaccumulation section using a 1D array.
   * \param[in] data - the output 1D array.
   */
  void SetPreaccOut(su2double* data, const int size);

  /*!
   * \brief Sets the input variables of a preaccumulation section using a 2D array.
   * \param[in] data - the output 1D array.
   */
  void SetPreaccOut(su2double** data, const int size_x, const int size_y);

  /*!
   * \brief Ends a preaccumulation section and computes the local Jacobi matrix
   * of a code section using the variables set with SetLocalInput(), SetLocalOutput() and pushes a statement
   * for each output variable to the AD tape.
   */
  void EndPreacc();

}

/*--- Macro to begin and end sections with a passive tape ---*/

#ifdef CODI_REVERSE_TYPE
#define AD_BEGIN_PASSIVE         \
  if(AD::globalTape.isActive()) {\
     AD::globalTape.setPassive();\
     AD::Status = true;          \
  }
#define AD_END_PASSIVE           \
  if(AD::Status) {               \
     AD::globalTape.setActive(); \
     AD::Status = false;         \
  }
#else
#define AD_BEGIN_PASSIVE
#define AD_END_PASSIVE
#endif

#include "ad_structure.inl"
