/*!
 * \file C1DInterpolation.hpp
 * \brief Classes for 1D interpolation.
 * \author Aman Baig, P. Gomes
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

#include <cassert>
#include <vector>
#include <algorithm>
#include "../option_structure.hpp"

/*!
 * \brief Base class for 1D interpolation.
 * \ingroup LookUpInterp
 */
class C1DInterpolation {
 protected:
  std::vector<su2double> x, y; /*!< \brief Data points. */

  /*!
   * \brief Find containing interval.
   * \param[in] xi - Point for which interpolation is required.
   * \return The start index of the interval, or the size if the coordinate is out of bounds.
   */
  inline size_t lower_bound(su2double xi) const {
    if (xi <= x.front() || xi >= x.back()) return x.size();

    size_t lb = 0, ub = x.size() - 1;

    while (ub - lb > 1) {
      size_t mid = (lb + ub) / 2;
      auto& change = (xi < x[mid]) ? ub : lb;
      change = mid;
    }
    return lb;
  }

 public:
  /*!
   * \brief Virtual destructor of the C1DInterpolation class.
   */
  virtual ~C1DInterpolation() = default;

  /*!
   * \brief Virtual method for setting the coefficients of the respective spline.
   * \param[in] X - the x values.
   * \param[in] Data - the f(x) values.
   */
  virtual void SetSpline(const std::vector<su2double>& X, const std::vector<su2double>& Data) {
    assert(X.size() == Data.size());
    x = X;
    y = Data;
  }

  /*!
   * \brief Evaluate the value of the spline at a point.
   */
  virtual su2double EvaluateSpline(su2double Point_Interp) const = 0;
  inline su2double operator()(su2double Point_Interp) const { return EvaluateSpline(Point_Interp); }
};

/*!
 * \brief Akima 1D interpolation.
 * \ingroup LookUpInterp
 */
class CAkimaInterpolation : public C1DInterpolation {
 protected:
  std::vector<su2double> b, c, d; /*!< \brief local variables for Akima spline cooefficients */

 public:
  CAkimaInterpolation() = default;

  /*!
   * \brief Constructor of the CAkimaInterpolation class.
   * \param[in] X - the x values (sorted low to high).
   * \param[in] Data - the f(x) values.
   */
  CAkimaInterpolation(const std::vector<su2double>& X, const std::vector<su2double>& Data) {
    CAkimaInterpolation::SetSpline(X, Data);
  }

  /*!
   * \brief Build the spline.
   */
  void SetSpline(const std::vector<su2double>& X, const std::vector<su2double>& Data) override;

  /*!
   * \brief Evaluate the value of the spline at a point.
   */
  su2double EvaluateSpline(su2double Point_Interp) const final;
};

/*!
 * \brief Cubic spline interpolation.
 * \ingroup LookUpInterp
 */
class CCubicSpline final : public CAkimaInterpolation {
 public:
  enum END_TYPE { SECOND, FIRST };

 private:
  const su2double startVal, endVal; /*!< \brief "boundary" values. */
  const END_TYPE startDer, endDer;  /*!< \brief 1st or 2nd derivative "boundary" conditions. */

 public:
  CCubicSpline() = default;

  /*!
   * \brief Constructor of the CCubicSpline class (defaults to natural spline).
   * \param[in] X - the x values (sorted low to high).
   * \param[in] Data - the f(x) values.
   * \param[in] startCondition - 1st or 2nd derivative imposed at the start.
   * \param[in] startValue - value of the derivative imposed at the start.
   * \param[in] endCondition - 1st or 2nd derivative imposed at the end.
   * \param[in] endValue - value of the derivative imposed at the end.
   */
  CCubicSpline(const std::vector<su2double>& X, const std::vector<su2double>& Data, END_TYPE startCondition = SECOND,
               su2double startValue = 0.0, END_TYPE endCondition = SECOND, su2double endValue = 0.0)
      : startVal(startValue), endVal(endValue), startDer(startCondition), endDer(endCondition) {
    SetSpline(X, Data);
  }

  /*!
   * \brief Build the spline.
   */
  void SetSpline(const std::vector<su2double>& X, const std::vector<su2double>& Data) override;
};

/*!
 * \brief Linear interpolation.
 * \ingroup LookUpInterp
 */
class CLinearInterpolation final : public C1DInterpolation {
 public:
  CLinearInterpolation() = default;

  /*!
   * \brief Constructor of the CLinearInterpolation class.
   * \param[in] X - the x values (sorted low to high).
   * \param[in] Data - the f(x) values.
   */
  CLinearInterpolation(const std::vector<su2double>& X, const std::vector<su2double>& Data) { SetSpline(X, Data); }

  /*!
   * \brief Evaluate the value of the spline at a point.
   */
  su2double EvaluateSpline(su2double Point_Interp) const override;
};

/*!
 * \brief Corrects for interpolation type.
 * \param[in] Inlet_Interpolated - the interpolated data after spline evaluation.
 * \param[in] Theta - the angle of the vertex (in xy plane).
 * \param[in] nDim - the dimensions of the case.
 * \param[in] Coord - the coordinates of the vertex.
 * \param[in] nVar_Turb - the number of turbulence variables as defined by turbulence model
 * \param[in] ENUM_INLET_INTERPOLATIONTYPE - enum of the interpolation type to be done
 * \returns the corrected Inlet Interpolated Data.
 */
std::vector<su2double> CorrectedInletValues(const std::vector<su2double>& Inlet_Interpolated, su2double Theta,
                                            unsigned short nDim, const su2double* Coord, unsigned short nVar_Turb,
                                            INLET_INTERP_TYPE Interpolation_Type);

/*!
 * \brief Prints the Inlet Interpolated Data
 * \param[in] Inlet_Interpolated_Interpolated - the final std::vector for the interpolated data
 * \param[in] Marker - name of the inlet marker
 * \param[in] nVertex - total number of vertexes.
 * \param[in] nDim - the dimensions of the problem.
 * \param[in] nColumns - the number of columns in the final interpolated data
 */
void PrintInletInterpolatedData(const std::vector<su2double>& Inlet_Data_Interpolated, const std::string& Marker,
                                unsigned long nVertex, unsigned short nDim, unsigned short nColumns);
