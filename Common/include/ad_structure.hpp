/*!
 * \file ad_structure.hpp
 * \brief Main routines for the algorithmic differentiation (AD) structure.
 * \author T. Albring
 * \version 6.2.0 "Falcon"
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
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
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
  void SetPreaccIn(const su2double* data, const int size);

  /*!
   * \brief Sets the input variables of a preaccumulation section using a 2D array.
   * \param[in] data - the input 2D array.
   * \param[in] size_x - size of the array in x dimension.
   * \param[in] size_y - size of the array in y dimension.
   */
  void SetPreaccIn(const su2double* const *data, const int size_x, const int size_y);

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
  
  /*!
   * \brief Initializes an externally differentiated function. Input and output variables are set with SetExtFuncIn/SetExtFuncOut
   * \param[in] storePrimalInput - Specifies whether the primal input values are stored for the reverse call of the external function.
   * \param[in] storePrimalOutput - Specifies whether the primal output values are stored for the reverse call of the external function.
   */
  void StartExtFunc(bool storePrimalInput, bool storePrimalOutput);
  
  /*!
   * \brief Sets the scalar input of a externally differentiated function.
   * \param[in] data - the scalar input variable.
   */
  void SetExtFuncIn(su2double &data);
  
  /*!
   * \brief Sets the input variables of a externally differentiated function using a 1D array.
   * \param[in] data - the input 1D array.
   * \param[in] size - number of rows.
   */
  void SetExtFuncIn(const su2double* data, const int size);
  
  /*!
  * \brief  Sets the input variables of a externally differentiated function using a 2D array.
  * \param[in] data - the input 2D array.
  * \param[in] size_x - number of rows.
  * \param[in] size_y - number of columns.
  */
  void SetExtFuncIn(const su2double* const *data, const int size_x, const int size_y);

  /*!
   * \brief Sets the scalar output of a externally differentiated function.
   * \param[in] data - the scalar output variable.
   */
  void SetExtFuncOut(su2double &data);

  /*!
   * \brief Sets the output variables of a externally differentiated function using a 1D array.
   * \param[in] data - the output 1D array.
   * \param[in] size - number of rows.
   */
  void SetExtFuncOut(su2double* data, const int size);

  /*!
  * \brief  Sets the output variables of a externally differentiated function using a 2D array.
  * \param[in] data - the output 2D array.
  * \param[in] size_x - number of rows.
  * \param[in] size_y - number of columns.
  */
  void SetExtFuncOut(su2double** data, const int size_x, const int size_y);
  
  /*!
   * \brief Ends an external function section by deleting the structures.
   */
  void EndExtFunc();

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

/*--- If we compile under OSX we have to overload some of the operators for
 *   complex numbers to avoid the use of the standard operators
 *  (they use a lot of functions that are only defined for doubles) ---*/

#ifdef __APPLE__

namespace std{
  template<>
  su2double abs(const complex<su2double>& x);

  template<>
  complex<su2double> operator/(const complex<su2double>& x, const complex<su2double>& y);

  template<>
  complex<su2double> operator*(const complex<su2double>& x, const complex<su2double>& y);
}
#endif
#include "ad_structure.inl"
