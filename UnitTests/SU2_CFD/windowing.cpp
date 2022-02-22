/*!
 * \file windowing.cpp
 * \brief Unit tests for windowed time-averaging.
 * \author C. Bauer
 * \version 7.3.0 "Blackbird"
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

#include "catch.hpp"
#include "../../SU2_CFD/include/output/COutput.hpp"
#include "../../SU2_CFD/include/output/tools/CWindowingTools.hpp"

/*!
 * \brief Wrapper class that provides reference values and methods for
 * averaging these over different windows.
 */
struct CWindowingTest {
  static constexpr unsigned long MAX_ITERATIONS = 1000;
  static constexpr unsigned long START_ITERATION = 10;
  static constexpr size_t MAX_INNER_ITERATIONS = 5;

  /*!
  * \brief Calculates sample value at the specified time-step. Simple SIN function.
  */
  static su2double GetSampleAtIteration(unsigned long currentTimeStep) {
    return 0.5 * sin(0.1 * currentTimeStep + 2) + 1.;
  }

  /*!
   * \brief Perform averaging over simulated inner and outer solver iterations.
   */
  static su2double calcAverage(WINDOW_FUNCTION win) {
    CWindowedAverage avg{win}; // Create averaging object with the specified window.
    su2double previousSample = 0.;
    su2double currentSample = 0.;
    su2double input = 0.;   // Value that is specified as actual input to the addValue function.
    su2double weightPrevious = 0.;
    su2double weightCurrent = 0.;
    // Simulate solver time-steps
    //printf("Checking window function %lu\n", (size_t)win);
    for (size_t outerIteration = 0; outerIteration < MAX_ITERATIONS; outerIteration++) {
      previousSample = currentSample;
      currentSample = GetSampleAtIteration(outerIteration);
      //printf("=== New timestep: %lu =======\n", outerIteration);
      //printf("PreviousSample: %f,\tCurrentSample: %f\n", previousSample, currentSample);
      // Simulate inner time-steps
      for (size_t innerIteration = 0; innerIteration < MAX_INNER_ITERATIONS; innerIteration++) {
        weightPrevious = static_cast<su2double>(MAX_INNER_ITERATIONS - innerIteration - 1);
        weightCurrent = innerIteration;
        // Simulate gradual change of sample during inner iterations.
        input = (weightPrevious * previousSample + weightCurrent * currentSample) / (weightCurrent + weightPrevious);
        avg.addValue(input, outerIteration, START_ITERATION);
        // printf("IT: %lu\twPrevious: %f\twCurrent: %f\tinput: %f\tavg: %f\n", 
        //   innerIteration, weightPrevious, weightCurrent, input, avg.GetVal());

      }
    }
    return avg.GetVal();
  }
  /*!
  * \brief Checks, if the specified result is within the expected limits
  */
  static bool checkResult(su2double result) { return (result < 1.1 && result > 0.9); }
};

TEST_CASE("BUMP", "[Windowing]") { 
    su2double avg = CWindowingTest::calcAverage(WINDOW_FUNCTION::BUMP);
    REQUIRE(CWindowingTest::checkResult(avg));
}

TEST_CASE("HANN", "[Windowing]") {
  su2double avg = CWindowingTest::calcAverage(WINDOW_FUNCTION::HANN);
  REQUIRE(CWindowingTest::checkResult(avg));
}

TEST_CASE("HANN_SQUARE", "[Windowing]") {
  su2double avg = CWindowingTest::calcAverage(WINDOW_FUNCTION::HANN_SQUARE);
  REQUIRE(CWindowingTest::checkResult(avg));
}

TEST_CASE("SQUARE", "[Windowing]") {
  su2double avg = CWindowingTest::calcAverage(WINDOW_FUNCTION::SQUARE);
  REQUIRE(CWindowingTest::checkResult(avg));
}
