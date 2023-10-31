/*!
 * \file CGradSmoothing.cpp
 * \brief Numerics for gradient smoothing problems.
 * \author T.Dick
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

#include "../../include/numerics/CGradSmoothing.hpp"

#include <limits>

CGradSmoothing::CGradSmoothing(unsigned short val_nDim, const CConfig* config) : CNumerics(val_nDim, val_nDim, config) {
  val_DHiDHj.resize(nDim, nDim);
  val_DHiDHj.setConstant(0.0);

  /*--- 8 is the max number of nodes in 3D ---*/
  Ni_Vec.resize(8);
  Ni_Vec.setConstant(0.0);
}

void CGradSmoothing::Compute_Tangent_Matrix(CElement *element, const CConfig *config) {

  unsigned int iDim, jDim, iGauss, nGauss, iShape, jShape, nNode;

  /*--- If we are on a curved design surface, everything is embedded in one dimension higher. --*/
  unsigned int nDimGlobal = nDim;
  if (config->GetSmoothOnSurface()) nDimGlobal=nDim+1;

  su2double Weight, Jac_X, val_HiHj, GradNiXGradNj = 0;

  su2double epsilon1 = config->GetSmoothingEps1();
  su2double epsilon2 = config->GetSmoothingEps2();

  /*--- Restarts the element: avoids adding over previous results in other elements --*/
  element->ClearElement();
  nNode = element->GetnNodes();
  nGauss = element->GetnGaussPoints();
  if (config->GetSmoothOnSurface()) {
    element->ComputeGrad_SurfaceEmbedded();
  } else {
    element->ComputeGrad_Linear();
  }

  /*--- Contribution from the gradients of the shape functions, representing the Laplace term. ---*/

  for (iGauss = 0; iGauss < nGauss; iGauss++) {

    Weight = element->GetWeight(iGauss);
    Jac_X = element->GetJ_X(iGauss);

    for (iShape = 0; iShape < nNode; iShape++) {

      /*--- Assumming symmetry ---*/
      for (jShape = iShape; jShape < nNode; jShape++) {

        for (iDim = 0; iDim < nDimGlobal; iDim++) {
          GradNiXGradNj += element->GetGradNi_X(iShape,iGauss,iDim)*element->GetGradNi_X(jShape,iGauss,iDim);
        }

        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0; jDim < nDim; jDim++) {
            if (iDim == jDim) {
              val_DHiDHj[iDim][jDim] = Weight * Jac_X * epsilon2 * GradNiXGradNj;
            } else {
              val_DHiDHj[iDim][jDim] = 0;
            }
          }
        }

        GradNiXGradNj=0;

        element->Add_DHiDHj(val_DHiDHj,iShape, jShape);
        /*--- Symmetric terms --*/
        if (iShape != jShape) {
          element->Add_DHiDHj_T(val_DHiDHj, jShape, iShape);
        }

      }

    }

  }

  /*--- Contribution from the shape functions themselves, representing the indentity term. --*/

  for (iGauss = 0; iGauss < nGauss; iGauss++) {

    Weight = element->GetWeight(iGauss);
    Jac_X = element->GetJ_X(iGauss);

    for (iShape = 0; iShape < nNode; iShape++) {
      Ni_Vec[iShape] = element->GetNi(iShape,iGauss);
    }

    for (iShape = 0; iShape < nNode; iShape++) {
      for (jShape = 0; jShape < nNode; jShape++) {
        val_HiHj = Weight * Jac_X * epsilon1 * Ni_Vec[iShape] * Ni_Vec[jShape];
        element->Add_HiHj(val_HiHj, iShape, jShape);
      }
    }

  }

}
