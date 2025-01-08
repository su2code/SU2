/*!
 * \file signal_processing_toolbox.hpp
 * \brief Header file for the signal processing toolbox.
 * \author S. Schotthöfer
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
  su2double cachedSum = 0.0;           /*!< \brief Cached sum of windowWeight*value over all previous iterations. */
  std::vector<su2double> values;       /*!< \brief Vector of instantatneous output values from starting time to the current time iteration.*/
  unsigned long lastTimeIter = std::numeric_limits<unsigned long>::max();
  const WINDOW_FUNCTION windowingFunctionId; /*!< \brief ID of the windowing function to use.*/

 public:

  /*!
   * \brief Creates a new CWindowedAverage with the specified windowing function
   */
  inline explicit CWindowedAverage(WINDOW_FUNCTION windowId) : windowingFunctionId(windowId)  {
    if (windowId==WINDOW_FUNCTION::SQUARE) {
      values.push_back(0.);
    }
  }

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
    cachedSum = 0.0;
    lastTimeIter = std::numeric_limits<unsigned long>::max();
  }

  /*!
   * \brief Adds the instantaneous output of the current iteration to the values-vector, if the current iteration is greater or equal to the starting iteration.
   * \param valIn - value of the instantaneous output, that should be added
   * \param currentIter - current time Iteration
   * \param startIter - iteration to start the windowed-time average.
   */
  void AddValue(su2double valIn, unsigned long curTimeIter,unsigned long startIter = 0);

private:
  /*!
  * \brief Caches the weighted sums from a previous time-step for later re-use
  * \param windowWidth - Total width of the window, over which the samples were weighted during the previous timestep
  */
  su2double UpdateCachedSum(unsigned long windowWidth) const;
};
