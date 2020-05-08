/*!
 * \file C1DInterpolation.hpp
 * \brief Inlet_interpolation_functions
 * \author Aman Baig
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include <iostream>
#include <cmath>
#include <vector>
#include<fstream>
#include "../datatype_structure.hpp"
#include "../option_structure.hpp"

using namespace std;

class C1DInterpolation{
protected:
  bool Point_Match = false; /*!< \brief to make sure all points from the inlet data match the vertex coordinates */
  vector<su2double> X;      /*!< \brief x for inlet data */
  vector<su2double> Data;   /*!< \brief f(x) for inlet data */

public:
  /*!
   * \brief Virtual destructor of the C1DInterpolation class.
   */
  virtual ~C1DInterpolation() = default;

  /*!
   * \brief virtual method for setting the cofficients for the respective spline spline.
   * \param[in] X - the x values.
   * \param[in] Data - the f(x) values.
   */
  virtual void SetSpline(const vector<su2double> &X, const vector<su2double> &Data){}

  /*!
   * \brief virtual method for evaluating the value of the respective Spline.
   * \param[in] Point_Interp - the point where interpolation value is required.
   * \returns the interpolated value.
   */
  virtual su2double EvaluateSpline(su2double Point_Interp){return 0;}

  /*!
   * \brief bool variable to make sure all vertex points fell in range of the inlet data.
   * \returns the bool variable.
   */
  bool GetPointMatch() const {return Point_Match;}
};

class CAkimaInterpolation final: public C1DInterpolation{
private:
  vector<su2double> x,y,b,c,d;  /*!< \brief local variables for Akima spline cooefficients */
  int n; /*!< \brief local variable for holding the size of the vector */

public:
  /*!
   * \brief Constructor of the CAkimaInterpolation class.
   * \param[in] X - the x values.
   * \param[in] Data - the f(x) values.
   */
  CAkimaInterpolation(vector<su2double> &X, vector<su2double> &Data){
      SetSpline(X,Data);
  }

  /*!
   * \brief for setting the cofficients for the Akima spline.
   * \param[in] X - the x values.
   * \param[in] Data - the f(x) values.
   */
  void SetSpline(const vector<su2double> &X, const vector<su2double> &Data) override;

  /*!
   * \brief For evaluating the value of Akima Spline.
   * \param[in] Point_Interp - the point where interpolation value is required.
   * \returns the interpolated value.
   */
  su2double EvaluateSpline(su2double Point_Interp) override;
};

class CLinearInterpolation final: public C1DInterpolation{
private:
  vector<su2double> x,y,dydx; /*!< \brief local variables for linear 'spline' cooefficients */
  int n;                      /*!< \brief local variable for holding the size of the vector */

public:
  /*!
   * \brief Constructor of the CLinearInterpolation class.
   * \param[in] X - the x values.
   * \param[in] Data - the f(x) values.
   */
  CLinearInterpolation(vector<su2double> &X, vector<su2double> &Data){
      SetSpline(X,Data);
  }

  /*!
   * \brief for setting the cofficients for Linear 'spline'.
   * \param[in] X - the x values.
   * \param[in] Data - the f(x) values.
   */
  void SetSpline(const vector<su2double> &X, const vector<su2double> &Data) override;

  /*!
   * \brief For evaluating the value for Linear 'spline'.
   * \param[in] Point_Interp - the point where interpolation value is required.
   * \returns the interpolated value.
   */
  su2double EvaluateSpline(su2double Point_Interp) override;
};

/*!
 * \brief to correct for interpolation type.
 * \param[in] Inlet_Interpolated - the interpolated data after spline evaluation.
 * \param[in] Theta - the angle of the vertex (in xy plane).
 * \param[in] nDim - the dimensions of the case.
 * \param[in] Coord - the coordinates of the vertex.
 * \param[in] nVar_Turb - the number of turbulence variables as defined by turbulence model
 * \param[in] ENUM_INLET_INTERPOLATIONTYPE - enum of the interpolation type to be done
 * \returns the corrected Inlet Interpolated Data.
 */
vector<su2double> CorrectedInletValues(const vector<su2double> &Inlet_Interpolated,
                                       su2double Theta ,
                                       unsigned short nDim,
                                       const su2double *Coord,
                                       unsigned short nVar_Turb,
                                       ENUM_INLET_INTERPOLATIONTYPE Interpolation_Type);

/*!
 * \brief to print the Inlet Interpolated Data
 * \param[in] Inlet_Interpolated_Interpolated - the final vector for the interpolated data
 * \param[in] Marker - name of the inlet marker
 * \param[in] nVertex - total number of vertexes.
 * \param[in] nDim - the dimensions of the problem.
 * \param[in] nColumns - the number of columns in the final interpolated data
 */
void PrintInletInterpolatedData(const vector<su2double>& Inlet_Data_Interpolated, string Marker,
                                unsigned long nVertex, unsigned short nDim, unsigned short nColumns);
