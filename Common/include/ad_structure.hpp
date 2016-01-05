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
   * \struct CArray1D
   * \brief The CArray1D struct is used to encapsulate a 1D array and its size
   *  so that it can be used as an argument for variadic functions.
   * \author T. Albring
   * \version 4.0.2 "Cardinal"
   */
  struct CArray1D{
    su2double* vec; /*!< \brief The actual pointer to the memory location. */
    unsigned short size; /*!< \brief The size of the array. */

    /*!
     * \brief Constructor of the struct.
     */
    CArray1D(su2double* vec_, unsigned short size_):
      vec(vec_), size(size_){}
  };

  /*!
   * \struct CArray2D
   * \brief The CArray2D struct is used to encapsulate a 2D array and its dimensions
   * so that it can be used as an argument for variadic functions.
   * \author T. Albring
   * \version 4.0.2 "Cardinal"
   */
  struct CArray2D{
    su2double** mat; /*!< \brief The actual pointer to the memory location. */
    unsigned short size_x, size_y; /*!< \brief The dimensions of the array. */

    /*!
     * \brief Constructor of the struct.
     */
    CArray2D(su2double** mat_, unsigned short size_x_, unsigned short size_y_):
      mat(mat_), size_x(size_x_), size_y(size_y_){}
  };

  /*!
   * \brief Sets the scalar input of a preaccumulation section.
   * \param[in] data - the scalar input variable.
   */
  void SetPreaccInput(const su2double &data);

  /*!
   * \brief Sets the input variables of a preaccumulation section using a 1D array.
   * \param[in] data - the input 1D array.
   */
  void SetPreaccInput(const CArray1D &data);

  /*!
   * \brief Sets the input variables of a preaccumulation section using a 2D array.
   * \param[in] data - the input 1D array.
   */
  void SetPreaccInput(const CArray2D &data);

  /*!
   * \brief Termination function for the variadic SetPreaccInput_Variadic() template.
   */
  void SetPreaccInput_Variadic();

  /*!
   * \brief Splits the argument list one by one into separate arguments
   *  and calls the corresponding SetPreaccInput() functions for su2double, CArray1D or CArray2D.
   * \param[in] arg1 - First argument.
   * \param[in] args - The remaining arguments.
   */
  template <typename Arg1, typename ... Args>
  void SetPreaccInput_Variadic(const Arg1& arg1, Args& ... args);

  /*!
   * \brief Starts a new preaccumulation section and sets the input variables. This function is variadic
   *  and accepts any combination of su2double, CArray1D and CArray2D.
   *
   * The idea of preaccumulation is to store only the Jacobi matrix of a code section during
   * the taping process instead of all operations. This decreases the tape size and reduces runtime.
   *
   * Typically a preaccumulation section has a lot of input variables. Instead of calling the same function
   * (SetPreaccInput()) for each variable we can pass all variables and arrays at once to this function.
   *
   * Example:
   * \code{.cpp}
   *  AD::StartPreacc(AD::CArray2D(TurbVar_Grad_i, nVar, nDim), AD::CArray1D(TurbVar_i, nVar), F1_i);
   * \endcode
   *
   * Note: the call of this routine must be followed by a call of EndPreacc() and the end of the code section.
   *
   * \param[in] args - The input variables of the preaccumulation section (any combination of su2double, CArray1D or CArray2D).
   *
   */
  template <typename ... Args>
  void StartPreacc(Args && ... args);

  /*!
   * \brief Sets the scalar output of a preaccumulation section.
   * \param[in] data - the scalar output variable.
   */
  void SetPreaccOutput(su2double &data);

  /*!
   * \brief Sets the output variables of a preaccumulation section using a 1D array.
   * \param[in] data - the output 1D array.
   */
  void SetPreaccOutput(CArray1D &data);

  /*!
   * \brief Sets the input variables of a preaccumulation section using a 2D array.
   * \param[in] data - the output 1D array.
   */
  void SetPreaccOutput(CArray2D &data);

  /*!
   * \brief Termination function for the variadic SetPreaccInput_Variadic() template.
   */
  void SetPreaccOutput_Variadic();

  /*!
   * \brief Splits the argument list one by one into separate arguments
   *  and calls the corresponding SetPreaccOutput() functions for su2double, CArray1D or CArray2D.
   * \param[in] arg1 - First argument.
   * \param[in] args - The remaining arguments.
   */
  template <typename Arg1, typename ... Args>
  void SetPreaccOutput_Variadic(Arg1& arg1, Args& ... args);

  /*!
   * \brief Ends a preaccumulation section and sets the output variables. This function is variadic
   *  and accepts any combination of su2double, CArray1D and CArray2D.
   *
   * \param[in] args - The output variables of the preaccumulation section (any combination of su2double, CArray1D or CArray2D).
   *
   */
  template <typename ... Args>
  void EndPreacc(Args && ... args);

  /*!
   * \brief Computes the local Jacobi matrix of a code section using the variables
   * set with SetLocalInput(), SetLocalOutput() and pushes a statement for each output variable
   * to the AD tape. Is called by EndPreacc().
   */
  void Preaccumulate();

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
