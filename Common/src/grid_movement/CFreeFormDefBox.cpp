/*!
 * \file CFreeFormDefBox.cpp
 * \brief Subroutines for handling Free-Form Deformation Boxes
 * \author F. Palacios, T. Economon, S. Padron
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

#include "../../include/grid_movement/CFreeFormDefBox.hpp"
#include "../../include/grid_movement/CBezierBlending.hpp"
#include "../../include/grid_movement/CBSplineBlending.hpp"
#include "../../include/toolboxes/geometry_toolbox.hpp"

CFreeFormDefBox::CFreeFormDefBox() : CGridMovement() {}

CFreeFormDefBox::CFreeFormDefBox(const unsigned short Degree[], unsigned short BSplineOrder[],
                                 unsigned short kind_blending)
    : CGridMovement() {
  unsigned short iCornerPoints, iOrder, jOrder, kOrder, iDim;

  /*--- FFD is always 3D (even in 2D problems) ---*/

  nDim = 3;
  nCornerPoints = 8;

  /*--- Allocate Corners points ---*/

  Coord_Corner_Points = new su2double*[nCornerPoints];
  for (iCornerPoints = 0; iCornerPoints < nCornerPoints; iCornerPoints++)
    Coord_Corner_Points[iCornerPoints] = new su2double[nDim];

  ParamCoord = new su2double[nDim];
  ParamCoord_ = new su2double[nDim];
  cart_coord = new su2double[nDim];
  cart_coord_ = new su2double[nDim];
  Gradient = new su2double[nDim];

  lDegree = Degree[0];
  lOrder = lDegree + 1;
  mDegree = Degree[1];
  mOrder = mDegree + 1;
  nDegree = Degree[2];
  nOrder = nDegree + 1;
  nControlPoints = lOrder * mOrder * nOrder;

  lDegree_Copy = Degree[0];
  lOrder_Copy = lDegree + 1;
  mDegree_Copy = Degree[1];
  mOrder_Copy = mDegree + 1;
  nDegree_Copy = Degree[2];
  nOrder_Copy = nDegree + 1;
  nControlPoints_Copy = lOrder_Copy * mOrder_Copy * nOrder_Copy;

  Coord_Control_Points = new su2double***[lOrder];
  ParCoord_Control_Points = new su2double***[lOrder];
  Coord_Control_Points_Copy = new su2double***[lOrder];
  for (iOrder = 0; iOrder < lOrder; iOrder++) {
    Coord_Control_Points[iOrder] = new su2double**[mOrder];
    ParCoord_Control_Points[iOrder] = new su2double**[mOrder];
    Coord_Control_Points_Copy[iOrder] = new su2double**[mOrder];
    for (jOrder = 0; jOrder < mOrder; jOrder++) {
      Coord_Control_Points[iOrder][jOrder] = new su2double*[nOrder];
      ParCoord_Control_Points[iOrder][jOrder] = new su2double*[nOrder];
      Coord_Control_Points_Copy[iOrder][jOrder] = new su2double*[nOrder];
      for (kOrder = 0; kOrder < nOrder; kOrder++) {
        Coord_Control_Points[iOrder][jOrder][kOrder] = new su2double[nDim];
        ParCoord_Control_Points[iOrder][jOrder][kOrder] = new su2double[nDim];
        Coord_Control_Points_Copy[iOrder][jOrder][kOrder] = new su2double[nDim];
        for (iDim = 0; iDim < nDim; iDim++) {
          Coord_Control_Points[iOrder][jOrder][kOrder][iDim] = 0.0;
          ParCoord_Control_Points[iOrder][jOrder][kOrder][iDim] = 0.0;
          Coord_Control_Points_Copy[iOrder][jOrder][kOrder][iDim] = 0.0;
        }
      }
    }
  }

  BlendingFunction = new CFreeFormBlending*[nDim];

  if (kind_blending == BEZIER) {
    BlendingFunction[0] = new CBezierBlending(lOrder, lOrder);
    BlendingFunction[1] = new CBezierBlending(mOrder, mOrder);
    BlendingFunction[2] = new CBezierBlending(nOrder, nOrder);
  }
  if (kind_blending == BSPLINE_UNIFORM) {
    BlendingFunction[0] = new CBSplineBlending(BSplineOrder[0], lOrder);
    BlendingFunction[1] = new CBSplineBlending(BSplineOrder[1], mOrder);
    BlendingFunction[2] = new CBSplineBlending(BSplineOrder[2], nOrder);
  }
}

CFreeFormDefBox::~CFreeFormDefBox() {
  unsigned short iOrder, jOrder, kOrder, iCornerPoints, iDim;

  for (iOrder = 0; iOrder < lOrder; iOrder++) {
    for (jOrder = 0; jOrder < mOrder; jOrder++) {
      for (kOrder = 0; kOrder < nOrder; kOrder++) {
        delete[] Coord_Control_Points[iOrder][jOrder][kOrder];
        delete[] ParCoord_Control_Points[iOrder][jOrder][kOrder];
        delete[] Coord_Control_Points_Copy[iOrder][jOrder][kOrder];
        if (Coord_SupportCP != nullptr) delete[] Coord_SupportCP[iOrder][jOrder][kOrder];
      }
      delete[] Coord_Control_Points[iOrder][jOrder];
      delete[] ParCoord_Control_Points[iOrder][jOrder];
      delete[] Coord_Control_Points_Copy[iOrder][jOrder];
      if (Coord_SupportCP != nullptr) delete[] Coord_SupportCP[iOrder][jOrder];
    }
    delete[] Coord_Control_Points[iOrder];
    delete[] ParCoord_Control_Points[iOrder];
    delete[] Coord_Control_Points_Copy[iOrder];
    if (Coord_SupportCP != nullptr) delete[] Coord_SupportCP[iOrder];
  }

  delete[] Coord_Control_Points;
  delete[] ParCoord_Control_Points;
  delete[] Coord_Control_Points_Copy;
  delete[] Coord_SupportCP;

  delete[] ParamCoord;
  delete[] ParamCoord_;
  delete[] cart_coord;
  delete[] cart_coord_;
  delete[] Gradient;

  for (iCornerPoints = 0; iCornerPoints < nCornerPoints; iCornerPoints++) delete[] Coord_Corner_Points[iCornerPoints];
  delete[] Coord_Corner_Points;

  for (iDim = 0; iDim < nDim; iDim++) {
    delete BlendingFunction[iDim];
  }
  delete[] BlendingFunction;
}

void CFreeFormDefBox::SetUnitCornerPoints() {
  unsigned short iDim;
  auto* coord = new su2double[nDim];

  for (iDim = 0; iDim < nDim; iDim++) coord[iDim] = 0.0;

  coord[0] = 0.0;
  coord[1] = 0.0;
  coord[2] = 0.0;
  this->SetCoordCornerPoints(coord, 0);
  coord[0] = 1.0;
  coord[1] = 0.0;
  coord[2] = 0.0;
  this->SetCoordCornerPoints(coord, 1);
  coord[0] = 1.0;
  coord[1] = 1.0;
  coord[2] = 0.0;
  this->SetCoordCornerPoints(coord, 2);
  coord[0] = 0.0;
  coord[1] = 1.0;
  coord[2] = 0.0;
  this->SetCoordCornerPoints(coord, 3);
  coord[0] = 0.0;
  coord[1] = 0.0;
  coord[2] = 1.0;
  this->SetCoordCornerPoints(coord, 4);
  coord[0] = 1.0;
  coord[1] = 0.0;
  coord[2] = 1.0;
  this->SetCoordCornerPoints(coord, 5);
  coord[0] = 1.0;
  coord[1] = 1.0;
  coord[2] = 1.0;
  this->SetCoordCornerPoints(coord, 6);
  coord[0] = 0.0;
  coord[1] = 1.0;
  coord[2] = 1.0;
  this->SetCoordCornerPoints(coord, 7);

  delete[] coord;
}

void CFreeFormDefBox::SetControlPoints_Parallelepiped() {
  unsigned short iDim, iDegree, jDegree, kDegree;

  /*--- Set base control points according to the notation of Vtk for hexahedrons ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    Coord_Control_Points[0][0][0][iDim] = Coord_Corner_Points[0][iDim];
    Coord_Control_Points[lOrder - 1][0][0][iDim] = Coord_Corner_Points[1][iDim];
    Coord_Control_Points[lOrder - 1][mOrder - 1][0][iDim] = Coord_Corner_Points[2][iDim];
    Coord_Control_Points[0][mOrder - 1][0][iDim] = Coord_Corner_Points[3][iDim];
    Coord_Control_Points[0][0][nOrder - 1][iDim] = Coord_Corner_Points[4][iDim];
    Coord_Control_Points[lOrder - 1][0][nOrder - 1][iDim] = Coord_Corner_Points[5][iDim];
    Coord_Control_Points[lOrder - 1][mOrder - 1][nOrder - 1][iDim] = Coord_Corner_Points[6][iDim];
    Coord_Control_Points[0][mOrder - 1][nOrder - 1][iDim] = Coord_Corner_Points[7][iDim];
  }

  /*--- Fill the rest of the cubic matrix of control points with uniform spacing (parallelepiped) ---*/
  for (iDegree = 0; iDegree <= lDegree; iDegree++)
    for (jDegree = 0; jDegree <= mDegree; jDegree++)
      for (kDegree = 0; kDegree <= nDegree; kDegree++) {
        Coord_Control_Points[iDegree][jDegree][kDegree][0] =
            Coord_Corner_Points[0][0] +
            su2double(iDegree) / su2double(lDegree) * (Coord_Corner_Points[1][0] - Coord_Corner_Points[0][0]);
        Coord_Control_Points[iDegree][jDegree][kDegree][1] =
            Coord_Corner_Points[0][1] +
            su2double(jDegree) / su2double(mDegree) * (Coord_Corner_Points[3][1] - Coord_Corner_Points[0][1]);
        Coord_Control_Points[iDegree][jDegree][kDegree][2] =
            Coord_Corner_Points[0][2] +
            su2double(kDegree) / su2double(nDegree) * (Coord_Corner_Points[4][2] - Coord_Corner_Points[0][2]);
      }
}

void CFreeFormDefBox::SetSupportCP(CFreeFormDefBox* FFDBox) {
  unsigned short iDim, iOrder, jOrder, kOrder;
  unsigned short lOrder = FFDBox->GetlOrder();
  unsigned short mOrder = FFDBox->GetmOrder();
  unsigned short nOrder = FFDBox->GetnOrder();

  Coord_SupportCP = new su2double***[lOrder];
  for (iOrder = 0; iOrder < lOrder; iOrder++) {
    Coord_SupportCP[iOrder] = new su2double**[mOrder];
    for (jOrder = 0; jOrder < mOrder; jOrder++) {
      Coord_SupportCP[iOrder][jOrder] = new su2double*[nOrder];
      for (kOrder = 0; kOrder < nOrder; kOrder++) Coord_SupportCP[iOrder][jOrder][kOrder] = new su2double[nDim];
    }
  }

  /*--- Set base support control points according to the notation of Vtk for hexahedrons ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    Coord_SupportCP[0][0][0][iDim] = Coord_Corner_Points[0][iDim];
    Coord_SupportCP[lOrder - 1][0][0][iDim] = Coord_Corner_Points[1][iDim];
    Coord_SupportCP[lOrder - 1][mOrder - 1][0][iDim] = Coord_Corner_Points[2][iDim];
    Coord_SupportCP[0][mOrder - 1][0][iDim] = Coord_Corner_Points[3][iDim];
    Coord_SupportCP[0][0][nOrder - 1][iDim] = Coord_Corner_Points[4][iDim];
    Coord_SupportCP[lOrder - 1][0][nOrder - 1][iDim] = Coord_Corner_Points[5][iDim];
    Coord_SupportCP[lOrder - 1][mOrder - 1][nOrder - 1][iDim] = Coord_Corner_Points[6][iDim];
    Coord_SupportCP[0][mOrder - 1][nOrder - 1][iDim] = Coord_Corner_Points[7][iDim];
  }

  /*--- Fill the rest of the cubic matrix of support control points with uniform spacing  ---*/
  for (iOrder = 0; iOrder < lOrder; iOrder++)
    for (jOrder = 0; jOrder < mOrder; jOrder++)
      for (kOrder = 0; kOrder < nOrder; kOrder++) {
        Coord_SupportCP[iOrder][jOrder][kOrder][0] =
            Coord_Corner_Points[0][0] +
            su2double(iOrder) / su2double(lOrder - 1) * (Coord_Corner_Points[1][0] - Coord_Corner_Points[0][0]);
        Coord_SupportCP[iOrder][jOrder][kOrder][1] =
            Coord_Corner_Points[0][1] +
            su2double(jOrder) / su2double(mOrder - 1) * (Coord_Corner_Points[3][1] - Coord_Corner_Points[0][1]);
        Coord_SupportCP[iOrder][jOrder][kOrder][2] =
            Coord_Corner_Points[0][2] +
            su2double(kOrder) / su2double(nOrder - 1) * (Coord_Corner_Points[4][2] - Coord_Corner_Points[0][2]);
      }
}

void CFreeFormDefBox::SetSupportCPChange(CFreeFormDefBox* FFDBox) {
  unsigned short iDim, iOrder, jOrder, kOrder;
  su2double *CartCoordNew, *ParamCoord;
  unsigned short lOrder = FFDBox->GetlOrder();
  unsigned short mOrder = FFDBox->GetmOrder();
  unsigned short nOrder = FFDBox->GetnOrder();

  auto**** ParamCoord_SupportCP = new su2double***[lOrder];
  for (iOrder = 0; iOrder < lOrder; iOrder++) {
    ParamCoord_SupportCP[iOrder] = new su2double**[mOrder];
    for (jOrder = 0; jOrder < mOrder; jOrder++) {
      ParamCoord_SupportCP[iOrder][jOrder] = new su2double*[nOrder];
      for (kOrder = 0; kOrder < nOrder; kOrder++) ParamCoord_SupportCP[iOrder][jOrder][kOrder] = new su2double[nDim];
    }
  }

  for (iOrder = 0; iOrder < lOrder; iOrder++)
    for (jOrder = 0; jOrder < mOrder; jOrder++)
      for (kOrder = 0; kOrder < nOrder; kOrder++)
        for (iDim = 0; iDim < nDim; iDim++)
          ParamCoord_SupportCP[iOrder][jOrder][kOrder][iDim] = Coord_SupportCP[iOrder][jOrder][kOrder][iDim];

  for (iDim = 0; iDim < nDim; iDim++) {
    Coord_Control_Points[0][0][0][iDim] = FFDBox->GetCoordCornerPoints(iDim, 0);
    Coord_Control_Points[1][0][0][iDim] = FFDBox->GetCoordCornerPoints(iDim, 1);
    Coord_Control_Points[1][1][0][iDim] = FFDBox->GetCoordCornerPoints(iDim, 2);
    Coord_Control_Points[0][1][0][iDim] = FFDBox->GetCoordCornerPoints(iDim, 3);
    Coord_Control_Points[0][0][1][iDim] = FFDBox->GetCoordCornerPoints(iDim, 4);
    Coord_Control_Points[1][0][1][iDim] = FFDBox->GetCoordCornerPoints(iDim, 5);
    Coord_Control_Points[1][1][1][iDim] = FFDBox->GetCoordCornerPoints(iDim, 6);
    Coord_Control_Points[0][1][1][iDim] = FFDBox->GetCoordCornerPoints(iDim, 7);
  }

  for (iOrder = 0; iOrder < FFDBox->GetlOrder(); iOrder++) {
    for (jOrder = 0; jOrder < FFDBox->GetmOrder(); jOrder++) {
      for (kOrder = 0; kOrder < FFDBox->GetnOrder(); kOrder++) {
        ParamCoord = ParamCoord_SupportCP[iOrder][jOrder][kOrder];
        CartCoordNew = EvalCartesianCoord(ParamCoord);
        FFDBox->SetCoordControlPoints(CartCoordNew, iOrder, jOrder, kOrder);
        FFDBox->SetCoordControlPoints_Copy(CartCoordNew, iOrder, jOrder, kOrder);
      }
    }
  }
}

void CFreeFormDefBox::SetCart2Cyl_ControlPoints(CConfig* config) {
  unsigned short iDegree, jDegree, kDegree;
  su2double CartCoord[3];
  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;

  X_0 = config->GetFFD_Axis(0);
  Y_0 = config->GetFFD_Axis(1);
  Z_0 = config->GetFFD_Axis(2);

  for (kDegree = 0; kDegree <= nDegree; kDegree++) {
    for (jDegree = 0; jDegree <= mDegree; jDegree++) {
      for (iDegree = 0; iDegree <= lDegree; iDegree++) {
        CartCoord[0] = Coord_Control_Points[iDegree][jDegree][kDegree][0];
        CartCoord[1] = Coord_Control_Points[iDegree][jDegree][kDegree][1];
        CartCoord[2] = Coord_Control_Points[iDegree][jDegree][kDegree][2];

        Xbar = CartCoord[0] - X_0;
        Ybar = CartCoord[1] - Y_0;
        Zbar = CartCoord[2] - Z_0;

        Coord_Control_Points[iDegree][jDegree][kDegree][0] = sqrt(Ybar * Ybar + Zbar * Zbar);
        Coord_Control_Points[iDegree][jDegree][kDegree][1] = atan2(Zbar, Ybar);
        if (Coord_Control_Points[iDegree][jDegree][kDegree][1] > PI_NUMBER / 2.0)
          Coord_Control_Points[iDegree][jDegree][kDegree][1] -= 2.0 * PI_NUMBER;
        Coord_Control_Points[iDegree][jDegree][kDegree][2] = Xbar;

        CartCoord[0] = Coord_Control_Points_Copy[iDegree][jDegree][kDegree][0];
        CartCoord[1] = Coord_Control_Points_Copy[iDegree][jDegree][kDegree][1];
        CartCoord[2] = Coord_Control_Points_Copy[iDegree][jDegree][kDegree][2];

        Xbar = CartCoord[0] - X_0;
        Ybar = CartCoord[1] - Y_0;
        Zbar = CartCoord[2] - Z_0;

        Coord_Control_Points_Copy[iDegree][jDegree][kDegree][0] = sqrt(Ybar * Ybar + Zbar * Zbar);
        Coord_Control_Points_Copy[iDegree][jDegree][kDegree][1] = atan2(Zbar, Ybar);
        if (Coord_Control_Points_Copy[iDegree][jDegree][kDegree][1] > PI_NUMBER / 2.0)
          Coord_Control_Points_Copy[iDegree][jDegree][kDegree][1] -= 2.0 * PI_NUMBER;
        Coord_Control_Points_Copy[iDegree][jDegree][kDegree][2] = Xbar;
      }
    }
  }
}

void CFreeFormDefBox::SetCyl2Cart_ControlPoints(CConfig* config) {
  unsigned short iDegree, jDegree, kDegree;
  su2double PolarCoord[3];

  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;
  X_0 = config->GetFFD_Axis(0);
  Y_0 = config->GetFFD_Axis(1);
  Z_0 = config->GetFFD_Axis(2);

  for (kDegree = 0; kDegree <= nDegree; kDegree++) {
    for (jDegree = 0; jDegree <= mDegree; jDegree++) {
      for (iDegree = 0; iDegree <= lDegree; iDegree++) {
        PolarCoord[0] = Coord_Control_Points[iDegree][jDegree][kDegree][0];
        PolarCoord[1] = Coord_Control_Points[iDegree][jDegree][kDegree][1];
        PolarCoord[2] = Coord_Control_Points[iDegree][jDegree][kDegree][2];

        Xbar = PolarCoord[2];
        Ybar = PolarCoord[0] * cos(PolarCoord[1]);
        Zbar = PolarCoord[0] * sin(PolarCoord[1]);

        PolarCoord[0] = Xbar + X_0;
        PolarCoord[1] = Ybar + Y_0;
        PolarCoord[2] = Zbar + Z_0;

        Coord_Control_Points[iDegree][jDegree][kDegree][0] = PolarCoord[0];
        Coord_Control_Points[iDegree][jDegree][kDegree][1] = PolarCoord[1];
        Coord_Control_Points[iDegree][jDegree][kDegree][2] = PolarCoord[2];
      }
    }
  }
}

void CFreeFormDefBox::SetCart2Cyl_CornerPoints(CConfig* config) {
  unsigned short iCornerPoint;
  su2double* CartCoord;
  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;

  X_0 = config->GetFFD_Axis(0);
  Y_0 = config->GetFFD_Axis(1);
  Z_0 = config->GetFFD_Axis(2);

  for (iCornerPoint = 0; iCornerPoint < 8; iCornerPoint++) {
    CartCoord = GetCoordCornerPoints(iCornerPoint);
    Xbar = CartCoord[0] - X_0;
    Ybar = CartCoord[1] - Y_0;
    Zbar = CartCoord[2] - Z_0;

    CartCoord[0] = sqrt(Ybar * Ybar + Zbar * Zbar);
    CartCoord[1] = atan2(Zbar, Ybar);
    if (CartCoord[1] > PI_NUMBER / 2.0) CartCoord[1] -= 2.0 * PI_NUMBER;
    CartCoord[2] = Xbar;
  }
}

void CFreeFormDefBox::SetCyl2Cart_CornerPoints(CConfig* config) {
  unsigned short iCornerPoint;
  su2double* PolarCoord;
  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;

  X_0 = config->GetFFD_Axis(0);
  Y_0 = config->GetFFD_Axis(1);
  Z_0 = config->GetFFD_Axis(2);

  for (iCornerPoint = 0; iCornerPoint < 8; iCornerPoint++) {
    PolarCoord = GetCoordCornerPoints(iCornerPoint);

    Xbar = PolarCoord[2];
    Ybar = PolarCoord[0] * cos(PolarCoord[1]);
    Zbar = PolarCoord[0] * sin(PolarCoord[1]);

    PolarCoord[0] = Xbar + X_0;
    PolarCoord[1] = Ybar + Y_0;
    PolarCoord[2] = Zbar + Z_0;
  }
}

void CFreeFormDefBox::SetCart2Sphe_ControlPoints(CConfig* config) {
  unsigned short iDegree, jDegree, kDegree;
  su2double CartCoord[3];
  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;

  X_0 = config->GetFFD_Axis(0);
  Y_0 = config->GetFFD_Axis(1);
  Z_0 = config->GetFFD_Axis(2);

  for (kDegree = 0; kDegree <= nDegree; kDegree++) {
    for (jDegree = 0; jDegree <= mDegree; jDegree++) {
      for (iDegree = 0; iDegree <= lDegree; iDegree++) {
        CartCoord[0] = Coord_Control_Points[iDegree][jDegree][kDegree][0];
        CartCoord[1] = Coord_Control_Points[iDegree][jDegree][kDegree][1];
        CartCoord[2] = Coord_Control_Points[iDegree][jDegree][kDegree][2];

        Xbar = CartCoord[0] - X_0;
        Ybar = CartCoord[1] - Y_0;
        Zbar = CartCoord[2] - Z_0;

        Coord_Control_Points[iDegree][jDegree][kDegree][0] = sqrt(Xbar * Xbar + Ybar * Ybar + Zbar * Zbar);
        Coord_Control_Points[iDegree][jDegree][kDegree][1] = atan2(Zbar, Ybar);
        if (Coord_Control_Points[iDegree][jDegree][kDegree][1] > PI_NUMBER / 2.0)
          Coord_Control_Points[iDegree][jDegree][kDegree][1] -= 2.0 * PI_NUMBER;
        Coord_Control_Points[iDegree][jDegree][kDegree][2] =
            acos(Xbar / Coord_Control_Points[iDegree][jDegree][kDegree][0]);

        CartCoord[0] = Coord_Control_Points_Copy[iDegree][jDegree][kDegree][0];
        CartCoord[1] = Coord_Control_Points_Copy[iDegree][jDegree][kDegree][1];
        CartCoord[2] = Coord_Control_Points_Copy[iDegree][jDegree][kDegree][2];

        Xbar = CartCoord[0] - X_0;
        Ybar = CartCoord[1] - Y_0;
        Zbar = CartCoord[2] - Z_0;

        Coord_Control_Points_Copy[iDegree][jDegree][kDegree][0] = sqrt(Xbar * Xbar + Ybar * Ybar + Zbar * Zbar);
        Coord_Control_Points_Copy[iDegree][jDegree][kDegree][1] = atan2(Zbar, Ybar);
        if (Coord_Control_Points_Copy[iDegree][jDegree][kDegree][1] > PI_NUMBER / 2.0)
          Coord_Control_Points_Copy[iDegree][jDegree][kDegree][1] -= 2.0 * PI_NUMBER;
        Coord_Control_Points_Copy[iDegree][jDegree][kDegree][2] =
            acos(Xbar / Coord_Control_Points_Copy[iDegree][jDegree][kDegree][0]);
      }
    }
  }
}

void CFreeFormDefBox::SetSphe2Cart_ControlPoints(CConfig* config) {
  unsigned short iDegree, jDegree, kDegree;
  su2double PolarCoord[3];

  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;
  X_0 = config->GetFFD_Axis(0);
  Y_0 = config->GetFFD_Axis(1);
  Z_0 = config->GetFFD_Axis(2);

  for (kDegree = 0; kDegree <= nDegree; kDegree++) {
    for (jDegree = 0; jDegree <= mDegree; jDegree++) {
      for (iDegree = 0; iDegree <= lDegree; iDegree++) {
        PolarCoord[0] = Coord_Control_Points[iDegree][jDegree][kDegree][0];
        PolarCoord[1] = Coord_Control_Points[iDegree][jDegree][kDegree][1];
        PolarCoord[2] = Coord_Control_Points[iDegree][jDegree][kDegree][2];

        Xbar = PolarCoord[0] * cos(PolarCoord[2]);
        Ybar = PolarCoord[0] * cos(PolarCoord[1]) * sin(PolarCoord[2]);
        Zbar = PolarCoord[0] * sin(PolarCoord[1]) * sin(PolarCoord[2]);

        PolarCoord[0] = Xbar + X_0;
        PolarCoord[1] = Ybar + Y_0;
        PolarCoord[2] = Zbar + Z_0;

        Coord_Control_Points[iDegree][jDegree][kDegree][0] = PolarCoord[0];
        Coord_Control_Points[iDegree][jDegree][kDegree][1] = PolarCoord[1];
        Coord_Control_Points[iDegree][jDegree][kDegree][2] = PolarCoord[2];
      }
    }
  }
}

void CFreeFormDefBox::SetCart2Sphe_CornerPoints(CConfig* config) {
  unsigned short iCornerPoint;
  su2double* CartCoord;
  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;

  X_0 = config->GetFFD_Axis(0);
  Y_0 = config->GetFFD_Axis(1);
  Z_0 = config->GetFFD_Axis(2);

  for (iCornerPoint = 0; iCornerPoint < 8; iCornerPoint++) {
    CartCoord = GetCoordCornerPoints(iCornerPoint);
    Xbar = CartCoord[0] - X_0;
    Ybar = CartCoord[1] - Y_0;
    Zbar = CartCoord[2] - Z_0;

    CartCoord[0] = sqrt(Xbar * Xbar + Ybar * Ybar + Zbar * Zbar);
    CartCoord[1] = atan2(Zbar, Ybar);
    if (CartCoord[1] > PI_NUMBER / 2.0) CartCoord[1] -= 2.0 * PI_NUMBER;
    CartCoord[2] = acos(Xbar / CartCoord[0]);
  }
}

void CFreeFormDefBox::SetSphe2Cart_CornerPoints(CConfig* config) {
  unsigned short iCornerPoint;
  su2double* PolarCoord;
  su2double X_0, Y_0, Z_0, Xbar, Ybar, Zbar;

  X_0 = config->GetFFD_Axis(0);
  Y_0 = config->GetFFD_Axis(1);
  Z_0 = config->GetFFD_Axis(2);

  for (iCornerPoint = 0; iCornerPoint < 8; iCornerPoint++) {
    PolarCoord = GetCoordCornerPoints(iCornerPoint);

    Xbar = PolarCoord[0] * cos(PolarCoord[2]);
    Ybar = PolarCoord[0] * cos(PolarCoord[1]) * sin(PolarCoord[2]);
    Zbar = PolarCoord[0] * sin(PolarCoord[1]) * sin(PolarCoord[2]);

    PolarCoord[0] = Xbar + X_0;
    PolarCoord[1] = Ybar + Y_0;
    PolarCoord[2] = Zbar + Z_0;
  }
}

void CFreeFormDefBox::SetCGNS(CGeometry* geometry, unsigned short iFFDBox, bool original) {
#ifdef HAVE_CGNS

  char FFDBox_filename[MAX_STRING_SIZE];
  bool new_file;
  unsigned short iDim, iDegree, jDegree, kDegree;
  unsigned short outDim;
  unsigned int pos;
  char zonename[33];
  int FFDBox_cgns_file;
  int cell_dim, phys_dim;
  int cgns_base = 0, cgns_family, cgns_zone, cgns_err, dummy;
  const char* basename;

  /*--- FFD output is always 3D (even in 2D problems),
   this is important for debuging ---*/
  nDim = geometry->GetnDim();
  cell_dim = nDim, phys_dim = 3;

  SPRINTF(FFDBox_filename, "ffd_boxes.cgns");

  new_file = (original) && (iFFDBox == 0);

  if (new_file) {
    cgns_err = cg_open(FFDBox_filename, CG_MODE_WRITE, &FFDBox_cgns_file);
    if (cgns_err) cg_error_print();
    cgns_err = cg_descriptor_write("Title", "Visualization of the FFD boxes generated by SU2_DEF.");
    if (cgns_err) cg_error_print();
  } else {
    cgns_err = cg_open(FFDBox_filename, CG_MODE_MODIFY, &FFDBox_cgns_file);
    if (cgns_err) cg_error_print();
  }

  if (original) {
    basename = "Original_FFD";
  } else {
    basename = "Deformed_FFD";
  }

  if (iFFDBox == 0) {
    cgns_err = cg_base_write(FFDBox_cgns_file, basename, cell_dim, phys_dim, &cgns_base);
    if (cgns_err) cg_error_print();
  }
  cgns_err = cg_family_write(FFDBox_cgns_file, cgns_base, Tag.c_str(), &cgns_family);
  if (cgns_err) cg_error_print();

  cgsize_t dims[9];
  dims[0] = lDegree + 1;
  dims[1] = mDegree + 1;
  if (cell_dim == 3) {
    dims[2] = nDegree + 1;
  }
  cgsize_t pointlen = 1;
  for (int ii = 0; ii < cell_dim; ii++) {
    dims[ii + cell_dim] = dims[ii] - 1;
    dims[ii + 2 * cell_dim] = 0;
    pointlen *= dims[ii];
  }

  auto* buffer = new passivedouble[pointlen];
  SPRINTF(zonename, "SU2_Zone_%d", SU2_TYPE::Int(iFFDBox));

  cgns_err = cg_zone_write(FFDBox_cgns_file, cgns_base, zonename, dims, CGNS_ENUMV(Structured), &cgns_zone);
  if (cgns_err) cg_error_print();
  cgns_err = cg_goto(FFDBox_cgns_file, cgns_base, zonename, 0, nullptr);
  if (cgns_err) cg_error_print();
  cgns_err = cg_famname_write(Tag.c_str());
  if (cgns_err) cg_error_print();

  const char* coord_names[3] = {"CoordinateX", "CoordinateY", "CoordinateZ"};
  for (iDim = 0; iDim < nDim; iDim++) {
    outDim = nDim == 2 ? 0 : nDegree;
    for (kDegree = 0; kDegree <= outDim; kDegree++) {
      for (jDegree = 0; jDegree <= mDegree; jDegree++) {
        for (iDegree = 0; iDegree <= lDegree; iDegree++) {
          pos = iDegree + jDegree * (lDegree + 1) + kDegree * (lDegree + 1) * (mDegree + 1);
          buffer[pos] = SU2_TYPE::GetValue(Coord_Control_Points[iDegree][jDegree][kDegree][iDim]);
        }
      }
    }
    cgns_err = cg_coord_write(FFDBox_cgns_file, cgns_base, cgns_zone, CGNS_ENUMV(RealDouble), coord_names[iDim], buffer,
                              &dummy);
    if (cgns_err) cg_error_print();
  }
  if (nDim == 2) {
    std::fill_n(buffer, pointlen, 0.0);
    cgns_err =
        cg_coord_write(FFDBox_cgns_file, cgns_base, cgns_zone, CGNS_ENUMV(RealDouble), coord_names[2], buffer, &dummy);
    if (cgns_err) cg_error_print();
  }

  delete[] buffer;
  cgns_err = cg_close(FFDBox_cgns_file);
  if (cgns_err) cg_error_print();
#else  // Not built with CGNS support
  cout << "CGNS file requested for FFD but SU2 was built without CGNS support. No file written"
       << "\n";
#endif
}

void CFreeFormDefBox::SetTecplot(CGeometry* geometry, unsigned short iFFDBox, bool original) {
  ofstream FFDBox_file;
  char FFDBox_filename[MAX_STRING_SIZE];
  bool new_file;
  unsigned short iDim, iDegree, jDegree, kDegree;

  /*--- FFD output is always 3D (even in 2D problems),
   this is important for debuging ---*/

  nDim = 3;

  SPRINTF(FFDBox_filename, "ffd_boxes.dat");

  new_file = (original) && (iFFDBox == 0);

  if (new_file) {
    FFDBox_file.open(FFDBox_filename, ios::out);
    FFDBox_file << "TITLE = \"Visualization of the FFD boxes generated by SU2_DEF.\"" << endl;
    if (nDim == 2)
      FFDBox_file << R"(VARIABLES = "x", "y")" << endl;
    else
      FFDBox_file << R"(VARIABLES = "x", "y", "z")" << endl;
  } else
    FFDBox_file.open(FFDBox_filename, ios::out | ios::app);

  FFDBox_file << "ZONE T= \"" << Tag;
  if (original)
    FFDBox_file << " (Original FFD)\"";
  else
    FFDBox_file << " (Deformed FFD)\"";
  if (nDim == 2)
    FFDBox_file << ", I=" << lDegree + 1 << ", J=" << mDegree + 1 << ", DATAPACKING=POINT" << endl;
  else
    FFDBox_file << ", I=" << lDegree + 1 << ", J=" << mDegree + 1 << ", K=" << nDegree + 1 << ", DATAPACKING=POINT"
                << endl;

  FFDBox_file.precision(15);

  if (nDim == 2) {
    for (jDegree = 0; jDegree <= mDegree; jDegree++) {
      for (iDegree = 0; iDegree <= lDegree; iDegree++) {
        for (iDim = 0; iDim < nDim; iDim++)
          FFDBox_file << scientific << Coord_Control_Points[iDegree][jDegree][0][iDim] << "\t";
        FFDBox_file << "\n";
      }
    }
  } else {
    for (kDegree = 0; kDegree <= nDegree; kDegree++) {
      for (jDegree = 0; jDegree <= mDegree; jDegree++) {
        for (iDegree = 0; iDegree <= lDegree; iDegree++) {
          for (iDim = 0; iDim < nDim; iDim++)
            FFDBox_file << scientific << Coord_Control_Points[iDegree][jDegree][kDegree][iDim] << "\t";
          FFDBox_file << "\n";
        }
      }
    }
  }

  FFDBox_file.close();
}

void CFreeFormDefBox::SetParaview(CGeometry* geometry, unsigned short iFFDBox, bool original) {
  ofstream FFDBox_file;
  char FFDBox_filename[MAX_STRING_SIZE];
  bool new_file;
  unsigned short iDim, iDegree, jDegree, kDegree;

  nDim = geometry->GetnDim();

  new_file = original;

  if (new_file)
    SPRINTF(FFDBox_filename, "ffd_boxes_%d.vtk", SU2_TYPE::Int(iFFDBox));
  else
    SPRINTF(FFDBox_filename, "ffd_boxes_def_%d.vtk", SU2_TYPE::Int(iFFDBox));

  FFDBox_file.open(FFDBox_filename, ios::out);
  FFDBox_file << "# vtk DataFile Version 3.0" << endl;
  FFDBox_file << "vtk output" << endl;
  FFDBox_file << "ASCII" << endl;
  FFDBox_file << "DATASET STRUCTURED_GRID" << endl;

  if (nDim == 2)
    FFDBox_file << "DIMENSIONS " << lDegree + 1 << " " << mDegree + 1 << " " << 1 << endl;
  else
    FFDBox_file << "DIMENSIONS " << lDegree + 1 << " " << mDegree + 1 << " " << nDegree + 1 << endl;
  if (nDim == 2)
    FFDBox_file << "POINTS " << (lDegree + 1) * (mDegree + 1) << " float" << endl;
  else
    FFDBox_file << "POINTS " << (lDegree + 1) * (mDegree + 1) * (nDegree + 1) << " float" << endl;

  FFDBox_file.precision(15);

  if (nDim == 2) {
    for (jDegree = 0; jDegree <= mDegree; jDegree++) {
      for (iDegree = 0; iDegree <= lDegree; iDegree++) {
        for (iDim = 0; iDim < nDim; iDim++)
          FFDBox_file << scientific << Coord_Control_Points[iDegree][jDegree][0][iDim] << "\t";
        FFDBox_file << " 0.0 \n";
      }
    }
  } else {
    for (kDegree = 0; kDegree <= nDegree; kDegree++) {
      for (jDegree = 0; jDegree <= mDegree; jDegree++) {
        for (iDegree = 0; iDegree <= lDegree; iDegree++) {
          for (iDim = 0; iDim < nDim; iDim++)
            FFDBox_file << scientific << Coord_Control_Points[iDegree][jDegree][kDegree][iDim] << "\t";
          FFDBox_file << "\n";
        }
      }
    }
  }

  FFDBox_file.close();
}

su2double* CFreeFormDefBox::GetParametricCoord_Analytical(const su2double* cart_coord) {
  unsigned short iDim;
  su2double *e1, *e2, *e3, *e12, *e23, *e13, *p;

  /*--- Auxiliary Basis Vectors of the deformed FFDBox ---*/
  e1 = new su2double[3];
  e2 = new su2double[3];
  e3 = new su2double[3];
  for (iDim = 0; iDim < nDim; iDim++) {
    e1[iDim] = Coord_Corner_Points[1][iDim] - Coord_Corner_Points[0][iDim];
    e2[iDim] = Coord_Corner_Points[3][iDim] - Coord_Corner_Points[0][iDim];
    e3[iDim] = Coord_Corner_Points[4][iDim] - Coord_Corner_Points[0][iDim];
  }

  /*--- Respective Cross-Products ---*/
  e12 = new su2double[3];
  e23 = new su2double[3];
  e13 = new su2double[3];
  CrossProduct(e1, e2, e12);
  CrossProduct(e1, e3, e13);
  CrossProduct(e2, e3, e23);

  /*--- p is Tranlated vector from the origin ---*/
  p = new su2double[3];
  for (iDim = 0; iDim < nDim; iDim++) p[iDim] = cart_coord[iDim] - Coord_Corner_Points[0][iDim];

  ParamCoord[0] = DotProduct(e23, p) / DotProduct(e23, e1);
  ParamCoord[1] = DotProduct(e13, p) / DotProduct(e13, e2);
  ParamCoord[2] = DotProduct(e12, p) / DotProduct(e12, e3);

  delete[] e1;
  delete[] e2;
  delete[] e3;
  delete[] e12;
  delete[] e23;
  delete[] e13;
  delete[] p;

  return ParamCoord;
}

su2double* CFreeFormDefBox::EvalCartesianCoord(su2double* ParamCoord) const {
  unsigned short iDim, iDegree, jDegree, kDegree;

  for (iDim = 0; iDim < nDim; iDim++) cart_coord[iDim] = 0.0;

  for (iDegree = 0; iDegree <= lDegree; iDegree++)
    for (jDegree = 0; jDegree <= mDegree; jDegree++)
      for (kDegree = 0; kDegree <= nDegree; kDegree++)
        for (iDim = 0; iDim < nDim; iDim++) {
          cart_coord[iDim] += Coord_Control_Points[iDegree][jDegree][kDegree][iDim] *
                              BlendingFunction[0]->GetBasis(iDegree, ParamCoord[0]) *
                              BlendingFunction[1]->GetBasis(jDegree, ParamCoord[1]) *
                              BlendingFunction[2]->GetBasis(kDegree, ParamCoord[2]);
        }

  return cart_coord;
}

su2double* CFreeFormDefBox::GetFFDGradient(su2double* val_coord, su2double* xyz) {
  unsigned short iDim, jDim, lmn[3];

  /*--- Set the Degree of the spline ---*/

  lmn[0] = lDegree;
  lmn[1] = mDegree;
  lmn[2] = nDegree;

  for (iDim = 0; iDim < nDim; iDim++) Gradient[iDim] = 0.0;

  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      Gradient[jDim] += GetDerivative2(val_coord, iDim, xyz, lmn) * GetDerivative3(val_coord, iDim, jDim, lmn);

  return Gradient;
}

void CFreeFormDefBox::GetFFDHessian(su2double* uvw, su2double* xyz, su2double** val_Hessian) {
  unsigned short iDim, jDim, lmn[3];

  /*--- Set the Degree of the spline ---*/

  lmn[0] = lDegree;
  lmn[1] = mDegree;
  lmn[2] = nDegree;

  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++) val_Hessian[iDim][jDim] = 0.0;

  /*--- Note that being all the functions linear combinations of polynomials, they are C^\infty,
   and the Hessian will be symmetric; no need to compute the under-diagonal part, for example ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    val_Hessian[0][0] += 2.0 * GetDerivative3(uvw, iDim, 0, lmn) * GetDerivative3(uvw, iDim, 0, lmn) +
                         GetDerivative2(uvw, iDim, xyz, lmn) * GetDerivative5(uvw, iDim, 0, 0, lmn);

    val_Hessian[1][1] += 2.0 * GetDerivative3(uvw, iDim, 1, lmn) * GetDerivative3(uvw, iDim, 1, lmn) +
                         GetDerivative2(uvw, iDim, xyz, lmn) * GetDerivative5(uvw, iDim, 1, 1, lmn);

    val_Hessian[2][2] += 2.0 * GetDerivative3(uvw, iDim, 2, lmn) * GetDerivative3(uvw, iDim, 2, lmn) +
                         GetDerivative2(uvw, iDim, xyz, lmn) * GetDerivative5(uvw, iDim, 2, 2, lmn);

    val_Hessian[0][1] += 2.0 * GetDerivative3(uvw, iDim, 0, lmn) * GetDerivative3(uvw, iDim, 1, lmn) +
                         GetDerivative2(uvw, iDim, xyz, lmn) * GetDerivative5(uvw, iDim, 0, 1, lmn);

    val_Hessian[0][2] += 2.0 * GetDerivative3(uvw, iDim, 0, lmn) * GetDerivative3(uvw, iDim, 2, lmn) +
                         GetDerivative2(uvw, iDim, xyz, lmn) * GetDerivative5(uvw, iDim, 0, 2, lmn);

    val_Hessian[1][2] += 2.0 * GetDerivative3(uvw, iDim, 1, lmn) * GetDerivative3(uvw, iDim, 2, lmn) +
                         GetDerivative2(uvw, iDim, xyz, lmn) * GetDerivative5(uvw, iDim, 1, 2, lmn);
  }

  val_Hessian[1][0] = val_Hessian[0][1];
  val_Hessian[2][0] = val_Hessian[0][2];
  val_Hessian[2][1] = val_Hessian[1][2];
}

su2double* CFreeFormDefBox::GetParametricCoord_Iterative(unsigned long iPoint, su2double* xyz,
                                                         const su2double* ParamCoordGuess, CConfig* config) {
  su2double *IndepTerm, SOR_Factor = 1.0, MinNormError, NormError, Determinant, AdjHessian[3][3],
                        Temp[3] = {0.0, 0.0, 0.0};
  unsigned short iDim, jDim, RandonCounter;
  unsigned long iter;

  su2double tol = config->GetFFD_Tol() * 1E-3;
  unsigned short it_max = config->GetnFFD_Iter();
  unsigned short Random_Trials = 500;

  /*--- Allocate the Hessian ---*/

  Hessian = new su2double*[nDim];
  IndepTerm = new su2double[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Hessian[iDim] = new su2double[nDim];
    ParamCoord[iDim] = ParamCoordGuess[iDim];
    IndepTerm[iDim] = 0.0;
  }

  RandonCounter = 0;
  MinNormError = 1E6;

  /*--- External iteration ---*/

  for (iter = 0; iter < (unsigned long)it_max * Random_Trials; iter++) {
    /*--- The independent term of the solution of our system is -Gradient(sol_old) ---*/

    Gradient = GetFFDGradient(ParamCoord, xyz);

    for (iDim = 0; iDim < nDim; iDim++) IndepTerm[iDim] = -Gradient[iDim];

    /*--- Hessian = The Matrix of our system, getHessian(sol_old,xyz,...) ---*/

    GetFFDHessian(ParamCoord, xyz, Hessian);

    /*--- Adjoint to Hessian ---*/

    AdjHessian[0][0] = Hessian[1][1] * Hessian[2][2] - Hessian[1][2] * Hessian[2][1];
    AdjHessian[0][1] = Hessian[0][2] * Hessian[2][1] - Hessian[0][1] * Hessian[2][2];
    AdjHessian[0][2] = Hessian[0][1] * Hessian[1][2] - Hessian[0][2] * Hessian[1][1];
    AdjHessian[1][0] = Hessian[1][2] * Hessian[2][0] - Hessian[1][0] * Hessian[2][2];
    AdjHessian[1][1] = Hessian[0][0] * Hessian[2][2] - Hessian[0][2] * Hessian[2][0];
    AdjHessian[1][2] = Hessian[0][2] * Hessian[1][0] - Hessian[0][0] * Hessian[1][2];
    AdjHessian[2][0] = Hessian[1][0] * Hessian[2][1] - Hessian[1][1] * Hessian[2][0];
    AdjHessian[2][1] = Hessian[0][1] * Hessian[2][0] - Hessian[0][0] * Hessian[2][1];
    AdjHessian[2][2] = Hessian[0][0] * Hessian[1][1] - Hessian[0][1] * Hessian[1][0];

    /*--- Determinant of Hessian ---*/

    Determinant =
        Hessian[0][0] * AdjHessian[0][0] + Hessian[0][1] * AdjHessian[1][0] + Hessian[0][2] * AdjHessian[2][0];

    /*--- Hessian inverse ---*/

    if (Determinant != 0) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Temp[iDim] = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          Temp[iDim] += AdjHessian[iDim][jDim] * IndepTerm[jDim] / Determinant;
        }
      }
      for (iDim = 0; iDim < nDim; iDim++) {
        IndepTerm[iDim] = Temp[iDim];
      }
    }

    /*--- Update with Successive over-relaxation ---*/

    for (iDim = 0; iDim < nDim; iDim++) {
      ParamCoord[iDim] = (1.0 - SOR_Factor) * ParamCoord[iDim] + SOR_Factor * (ParamCoord[iDim] + IndepTerm[iDim]);
    }

    /*--- If the gradient is small, we have converged ---*/

    if ((fabs(IndepTerm[0]) < tol) && (fabs(IndepTerm[1]) < tol) && (fabs(IndepTerm[2]) < tol)) break;

    /*--- Compute the norm of the error ---*/

    NormError = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) NormError += IndepTerm[iDim] * IndepTerm[iDim];
    NormError = sqrt(NormError);

    MinNormError = min(NormError, MinNormError);

    /*--- If we have no convergence with Random_Trials iterations probably we are in a local minima. ---*/

    if (((iter % it_max) == 0) && (iter != 0)) {
      RandonCounter++;
      if (RandonCounter == Random_Trials) {
        cout << endl
             << "Unknown point: " << iPoint << " (" << xyz[0] << ", " << xyz[1] << ", " << xyz[2]
             << "). Min Error: " << MinNormError << ". Iter: " << iter << "." << endl;
      } else {
        SOR_Factor = 0.1;
        for (iDim = 0; iDim < nDim; iDim++) ParamCoord[iDim] = su2double(rand()) / su2double(RAND_MAX);
      }
    }

    /* --- Splines are not defined outside of [0,1]. So if the parametric coords are outside of
     *  [0,1] the step was too big and we have to use a smaller relaxation factor. ---*/

    if ((config->GetFFD_Blending() == BSPLINE_UNIFORM) &&
        (((ParamCoord[0] < 0.0) || (ParamCoord[0] > 1.0)) || ((ParamCoord[1] < 0.0) || (ParamCoord[1] > 1.0)) ||
         ((ParamCoord[2] < 0.0) || (ParamCoord[2] > 1.0)))) {
      for (iDim = 0; iDim < nDim; iDim++) {
        ParamCoord[iDim] = ParamCoordGuess[iDim];
      }
      SOR_Factor = 0.9 * SOR_Factor;
    }
  }

  for (iDim = 0; iDim < nDim; iDim++) delete[] Hessian[iDim];
  delete[] Hessian;
  delete[] IndepTerm;

  /*--- The code has hit the max number of iterations ---*/

  if (iter == (unsigned long)it_max * Random_Trials) {
    cout << "Unknown point: (" << xyz[0] << ", " << xyz[1] << ", " << xyz[2]
         << "). Increase the value of FFD_ITERATIONS." << endl;
  }

  /*--- Real Solution is now ParamCoord; Return it ---*/

  return ParamCoord;
}

bool CFreeFormDefBox::CheckPointInsideFFD(const su2double* coord) const {
  /*--- Indices of the FFD box. Note that the front face is labelled 0,1,2,3 and the back face is 4,5,6,7 ---*/

  unsigned short Index[6][5] = {{0, 1, 2, 3, 0},   // front side
                                {1, 5, 6, 2, 1},   // right side
                                {2, 6, 7, 3, 2},   // top side
                                {3, 7, 4, 0, 3},   // left side
                                {4, 5, 1, 0, 4},   // bottom side
                                {4, 7, 6, 5, 4}};  // back side

  /*--- The current approach is to subdivide each of the 6 faces of the hexahedral FFD box into 4 triangles by defining
   a supporting middle point. This allows nonplanar FFD boxes. Note that the definition of the FFD box is as follows:
   the FFD box is a 6-sided die and we are looking at the side "1". The opposite side is side "6". If we are looking at
   side "1", we define the nodes counterclockwise. If we are looking at side "6", we define the face clockwise. ---*/

  /*--- Loop over the faces of the FFD box. ---*/

  for (unsigned short iFace = 0; iFace < 6; iFace++) {
    /*--- Every face needs an interpolated middle point for the triangles. ---*/

    su2double P[3] = {0.0, 0.0, 0.0};
    for (int p = 0; p < 4; p++) {
      P[0] += 0.25 * Coord_Corner_Points[Index[iFace][p]][0];
      P[1] += 0.25 * Coord_Corner_Points[Index[iFace][p]][1];
      P[2] += 0.25 * Coord_Corner_Points[Index[iFace][p]][2];
    }

    /*--- Loop over the 4 triangles making up the FFD box. The sign should be equal for all distances. ---*/

    for (int iNode = 0; iNode < 4; iNode++) {
      const su2double* plane[] = {P, Coord_Corner_Points[Index[iFace][iNode]],
                                  Coord_Corner_Points[Index[iFace][iNode + 1]]};
      if (GeometryToolbox::PointToPlaneDistance(plane, coord) < 0) {
        return false;
      }
    }
  }
  return true;
}

su2double CFreeFormDefBox::GetDerivative1(su2double* uvw, unsigned short val_diff, unsigned short* ijk,
                                          unsigned short* lmn) const {
  unsigned short iDim;
  su2double value = 0.0;

  value = BlendingFunction[val_diff]->GetDerivative(ijk[val_diff], uvw[val_diff], 1);
  for (iDim = 0; iDim < nDim; iDim++)
    if (iDim != val_diff) value *= BlendingFunction[iDim]->GetBasis(ijk[iDim], uvw[iDim]);

  return value;
}

su2double CFreeFormDefBox::GetDerivative2(su2double* uvw, unsigned short dim, const su2double* xyz,
                                          const unsigned short* lmn) const {
  unsigned short iDegree, jDegree, kDegree;
  su2double value = 0.0;

  for (iDegree = 0; iDegree <= lmn[0]; iDegree++)
    for (jDegree = 0; jDegree <= lmn[1]; jDegree++)
      for (kDegree = 0; kDegree <= lmn[2]; kDegree++) {
        value += Coord_Control_Points[iDegree][jDegree][kDegree][dim] * BlendingFunction[0]->GetBasis(iDegree, uvw[0]) *
                 BlendingFunction[1]->GetBasis(jDegree, uvw[1]) * BlendingFunction[2]->GetBasis(kDegree, uvw[2]);
      }

  return 2.0 * (value - xyz[dim]);
}

su2double CFreeFormDefBox::GetDerivative3(su2double* uvw, unsigned short dim, unsigned short diff_this,
                                          unsigned short* lmn) const {
  unsigned short iDegree, jDegree, kDegree, iDim;
  su2double value = 0;

  auto* ijk = new unsigned short[nDim];

  for (iDim = 0; iDim < nDim; iDim++) ijk[iDim] = 0;

  for (iDegree = 0; iDegree <= lmn[0]; iDegree++)
    for (jDegree = 0; jDegree <= lmn[1]; jDegree++)
      for (kDegree = 0; kDegree <= lmn[2]; kDegree++) {
        ijk[0] = iDegree;
        ijk[1] = jDegree;
        ijk[2] = kDegree;
        value += Coord_Control_Points[iDegree][jDegree][kDegree][dim] * GetDerivative1(uvw, diff_this, ijk, lmn);
      }

  delete[] ijk;

  return value;
}

su2double CFreeFormDefBox::GetDerivative4(su2double* uvw, unsigned short val_diff, unsigned short val_diff2,
                                          unsigned short* ijk, unsigned short* lmn) const {
  unsigned short iDim;
  su2double value = 0.0;

  if (val_diff == val_diff2) {
    value = BlendingFunction[val_diff]->GetDerivative(ijk[val_diff], uvw[val_diff], 2);
    for (iDim = 0; iDim < nDim; iDim++)
      if (iDim != val_diff) value *= BlendingFunction[iDim]->GetBasis(ijk[iDim], uvw[iDim]);
  } else {
    value = BlendingFunction[val_diff]->GetDerivative(ijk[val_diff], uvw[val_diff], 1) *
            BlendingFunction[val_diff2]->GetDerivative(ijk[val_diff2], uvw[val_diff2], 1);
    for (iDim = 0; iDim < nDim; iDim++)
      if ((iDim != val_diff) && (iDim != val_diff2)) value *= BlendingFunction[iDim]->GetBasis(ijk[iDim], uvw[iDim]);
  }

  return value;
}

su2double CFreeFormDefBox::GetDerivative5(su2double* uvw, unsigned short dim, unsigned short diff_this,
                                          unsigned short diff_this_also, unsigned short* lmn) const {
  unsigned short iDegree, jDegree, kDegree, iDim;
  su2double value = 0.0;

  auto* ijk = new unsigned short[nDim];

  for (iDim = 0; iDim < nDim; iDim++) ijk[iDim] = 0;

  for (iDegree = 0; iDegree <= lmn[0]; iDegree++)
    for (jDegree = 0; jDegree <= lmn[1]; jDegree++)
      for (kDegree = 0; kDegree <= lmn[2]; kDegree++) {
        ijk[0] = iDegree;
        ijk[1] = jDegree;
        ijk[2] = kDegree;
        value += Coord_Control_Points[iDegree][jDegree][kDegree][dim] *
                 GetDerivative4(uvw, diff_this, diff_this_also, ijk, lmn);
      }

  delete[] ijk;

  return value;
}
