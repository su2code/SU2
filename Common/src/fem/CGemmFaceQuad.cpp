/*!
 * \file CGemmFaceQuad.cpp
 * \brief Functions for the class CGemmFaceQuad.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
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

#include "../../include/fem/CGemmFaceQuad.hpp"

/*----------------------------------------------------------------------------------*/
/*                  Public member functions of CGemmFaceQuad.                       */
/*----------------------------------------------------------------------------------*/

CGemmFaceQuad::CGemmFaceQuad(const int val_M, const int val_Type, const int val_K)
  : CGemmBase() {

  /*--- Copy the arguments into the member variables. ---*/
  M = val_M;
  K = val_K;

  TypeTensorProduct = val_Type;

  /*--- Determine the type of tensor product to be carried out. ---*/
  switch( TypeTensorProduct ) {

    case CGemmBase::DOFS_TO_INT: {

      /*--- Tensor product to create the data in the integration points of the
            line from the DOFs of the adjacent quadrilateral. Create the map
            with function pointers for this tensor product. ---*/
      map<CUnsignedShort2T, TPIS2D> mapFunctions;
      CreateMapTensorProductSurfaceIntPoints2D(mapFunctions);

      /*--- Try to find the combination of K and M in mapFunctions. If not found,
            write a clear error message that this tensor product is not supported. ---*/
      CUnsignedShort2T MK(M, K);
      auto MI = mapFunctions.find(MK);
      if(MI == mapFunctions.end()) {
        std::ostringstream message;
        message << "The tensor product TensorProductSurfaceIntPoints2D_" << M
                << "_" << K << " not created by the automatic source code "
                << "generator. Modify this automatic source code creator";
        SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
      }

      /*--- Set the function pointer to carry out tensor product. ---*/
      TensorProductDataSurfIntPoints = MI->second;

      break;
    }

    case CGemmBase::INT_TO_DOFS: {

      /*--- Tensor product to update the residuals of the adjacent quadrilateral
            from the data in the integration points of the line. Create the map
            with function pointers for this tensor product. ---*/
      map<CUnsignedShort2T, TPDR2D> mapFunctions;
      CreateMapTensorProductSurfaceResVolumeDOFs2D(mapFunctions);

      /*--- Try to find the combination of K and M in mapFunctions. If not found,
            write a clear error message that this tensor product is not supported. ---*/
      CUnsignedShort2T KM(K, M);
      auto MI = mapFunctions.find(KM);
      if(MI == mapFunctions.end()) {
        std::ostringstream message;
        message << "The tensor product TensorProductSurfaceResVolumeDOFs2D_" << K
                << "_" << M << " not created by the automatic source code "
                << "generator. Modify this automatic source code creator";
        SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
      }

      /*--- Set the function pointer to carry out tensor product. ---*/
      TensorProductResSurfIntPoints = MI->second;

      break;
    }

    default:
      SU2_MPI::Error(string("Invalid value of TypeTensorProduct"), CURRENT_FUNCTION);

  }
}

void CGemmFaceQuad::DOFs2Int(vector<ColMajorMatrix<passivedouble> > &tensor,
                             const int                              faceID_Elem,
                             const int                              N,
                             ColMajorMatrix<su2double>              &dataDOFs,
                             ColMajorMatrix<su2double>              &dataInt) {

  TensorProductDataSurfIntPoints(N, faceID_Elem, dataDOFs.rows(), dataInt.rows(),
                                 tensor[0].data(), tensor[1].data(), dataDOFs.data(),
                                 dataInt.data());
}

void CGemmFaceQuad::Int2DOFs(vector<ColMajorMatrix<passivedouble> > &tensor,
                             const int                              faceID_Elem,
                             const int                              N,
                             ColMajorMatrix<su2double>              &dataInt,
                             ColMajorMatrix<su2double>              &dataDOFs) {

  TensorProductResSurfIntPoints(N, faceID_Elem, dataInt.rows(), dataDOFs.rows(),
                                tensor[0].data(), tensor[1].data(), dataInt.data(),
                                dataDOFs.data());
}
