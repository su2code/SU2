/*!
 * \file signal_processing_toolbox.hpp
 * \brief Header file for the signal processing toolbox.
 * \author S. Schotthöfer
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include <vector>
#include <limits>
#include "../../../../Common/include/option_structure.hpp"

class CWindowingTools{
public:
  /*! \brief Returns the value of a windowing function given by windowId at time-step curTimeIter with given time-frame endTimeIter (i.e. endTimeIter=nTimeIter, if one starts  windowing at time t =0.)
   * \param windowId - enum specifing the used window
   * \param curTimeIter - Current time iteration of the solver
   * \param endTimeIter - Number of time steps to average over
   * \return Value of the window-weight-function at time curTimeIter with time-frame endTimeIter
   */
  static su2double GetWndWeight(WINDOW_FUNCTION windowId, unsigned long curTimeIter, unsigned long endTimeIter);

protected:
  // Long time windows
  /*! \brief Returns the value of the Hann-window function at time-step curTimeIter with given end-time endTimeIter.
  * \param curTimeIter - Current time iteration of the solver
  * \param endTimeIter - Number of time steps to average over
  * \return Value of the window-weight-function at time curTimeIter with end-time endTimeIter
  */
  static su2double HannWindow(unsigned long curTimeIter, unsigned long endTimeIter);

  /*! \brief Returns the value of the Hann-Square-window function at time-step i with given end-time endTimeIter.
  * \param curTimeIter - Current time iteration of the solver
  * \param endTimeIter - Number of time steps to average over
  * \return Value of the window-weight-function at time curTimeIter with end-time endTimeIter
  */
  static su2double HannSquaredWindow(unsigned long curTimeIter, unsigned long endTimeIter);

  /*! \brief Returns the value of the Bump-window function at time-step i with given end-time endTimeIter.
  * \param curTimeIter - Current time iteration of the solver
  * \param endTimeIter - Number of time steps to average over
  * \return Value of the window-weight-function at time curTimeIter with end-time endTimeIter
  */
  static su2double BumpWindow(unsigned long curTimeIter, unsigned long endTimeIter);
};

class CWindowedAverage:CWindowingTools{
private:
  su2double val = 0.0;                 /*!< \brief Value of the windowed-time average (of the instantaneous output) from starting time to the current time iteration. */
  std::vector<su2double> values;       /*!< \brief Vector of instantatneous output values from starting time to the current time iteration.*/
  unsigned long lastTimeIter = std::numeric_limits<unsigned long>::max();

public:
  /*!
   * \brief Returns the value of windowed-time average (of the instantaneous output) from starting time to the current time iteration
   */
  inline su2double GetVal() const { return val; }

  /*!
   *\brief Resets the value of windowed-time average (of the instantaneous output) from starting time to the current time iteration to 0.
   */
  inline void Reset() {
    val = 0.0;
    values.clear();
    lastTimeIter = std::numeric_limits<unsigned long>::max();
  }

  /*!
   * \brief Adds the instantaneous output of the current iteration to the values-vector, if the current iteration is greater or equal to the starting iteration.
   * \param valIn - value of the instantaneous output, that should be added
   * \param currentIter - current time Iteration
   * \param startIter - iteration to start the windowed-time average.
   */
  void addValue(su2double valIn, unsigned long curTimeIter,unsigned long startIter = 0);

  /*!
   * \brief Computes a windowed-time average of the values stored in the vector "values" using the windowing-function specified in enum windowId
   *        and stores it in "val".
   * \param windowId - specified windowing-function
   * \return windowed-time average of the values stored in the vector "values"
   */
  su2double WindowedUpdate(WINDOW_FUNCTION windowId);

private:
  /*! \brief Computes a Square-windowed-time average of the values stored in the vector "values" with the Midpoint-integration rule (for consistency with the adjoint solver).
  * \return  Squarewindowed-time average of the values stored in the vector "values"
  */
  su2double NoWindowing() const;

  /*! \brief Computes a Hann-windowed-time average of the values stored in the vector "values" with the Midpoint-integration rule (for consistency with the adjoint solver).
  * \return  Squarewindowed-time average of the values stored in the vector "values"
  */
  su2double HannWindowing() const;

  /*! \brief Computes a Hann-Square-windowed-time average of the values stored in the vector "values" with the Midpoint-integration rule (for consistency with the adjoint solver).
  * \return  Squarewindowed-time average of the values stored in the vector "values"
  */
  su2double HannSquaredWindowing() const;

  /*! \brief Computes a Bump-windowed-time average of the values stored in the vector "values" with the Midpoint-integration rule (for consistency with the adjoint solver).
  * \return  Squarewindowed-time average of the values stored in the vector "values"
  */
  su2double BumpWindowing() const;
};
