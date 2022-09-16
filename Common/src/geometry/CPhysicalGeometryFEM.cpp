/*!
 * \file CPhysicalGeometry.cpp
 * \brief Implementation of the FEM part physical geometry class.
 * \author F. Palacios, T. Economon
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

#include "../../include/toolboxes/CLinearPartitioner.hpp"
#include "../../include/toolboxes/fem/CMatchingFace.hpp"
#include "../../include/geometry/primal_grid/CPrimalGridFEM.hpp"
#include "../../include/geometry/primal_grid/CPrimalGridBoundFEM.hpp"
#include "../../include/geometry/CPhysicalGeometry.hpp"
#include "../../include/fem/CFEMStandardTriPartition.hpp"
#include "../../include/fem/CFEMStandardQuadPartition.hpp"
#include "../../include/fem/CFEMStandardTetPartition.hpp"
#include "../../include/fem/CFEMStandardPyraPartition.hpp"
#include "../../include/fem/CFEMStandardPrismPartition.hpp"
#include "../../include/fem/CFEMStandardHexPartition.hpp"
#include "../../include/fem/CFEMStandardLinePartition.hpp"

void CPhysicalGeometry::LoadLinearlyPartitionedPointsFEM(CConfig *config, CMeshReaderFVM *mesh) {

  /*--- Get the partitioned coordinates and their global IDs from the mesh object. ---*/
  const auto &gridCoords     = mesh->GetLocalPointCoordinates();
  const auto &globalPointIDs = mesh->GetGlobalPointIDs();

  /*--- Initialize point counts and the grid node data structure. ---*/

  nodes = new CPoint(nPoint, nDim);

  /*--- Loop over the points and set the coordinates and global index. ---*/
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    for (unsigned short iDim = 0; iDim < nDim; ++iDim)
      nodes->SetCoord(iPoint, iDim, gridCoords[iDim][iPoint]);
    nodes->SetGlobalIndex(iPoint, globalPointIDs[iPoint]);
  }
}

void CPhysicalGeometry::LoadLinearlyPartitionedVolumeElementsFEM(CConfig *config, CMeshReaderFVM *mesh) {

  /*--- Reset the global to local element mapping. ---*/
  Global_to_Local_Elem.clear();

  /*--- Get the volume connectivity from the mesh object. ---*/
  const auto &dataElems = mesh->GetLocalVolumeElementConnectivity();

  /*--- Allocate space for the interior elements in our SU2 data
        structure. Note that we only instantiate our rank's local set. ---*/
  elem = new CPrimalGrid*[nElem] ();

  /*--- Loop over all of the internal, local volumetric elements. ---*/
  unsigned long ind = 0;
  for (unsigned long jElem=0; jElem<nElem; ++jElem) {

    /*--- Create a FEM element from the data dataElems. ---*/
    const auto *dataElem = dataElems.data() + ind;
    elem[jElem] = new CPrimalGridFEM(dataElem);

    /*--- Store the global to local mapping in Global_to_Local_Elem. ---*/
    Global_to_Local_Elem[dataElem[4]] = jElem;

    /*--- Update ind for the next element. ---*/
    ind += dataElem[3] + 5;
  }
}

void CPhysicalGeometry::LoadLinearlyPartitionedSurfaceElementsFEM(CConfig *config, CMeshReaderFVM *mesh) {

  /*--- Store the number of markers and print to the screen. ---*/
  nMarker = mesh->GetNumberOfMarkers();
  config->SetnMarker_All(nMarker);
  if (rank == MASTER_NODE)
    cout << nMarker << " surface markers." << endl;

  /*--- Create the data structure for boundary elements. ---*/
  bound         = new CPrimalGrid**[nMarker];
  nElem_Bound   = new unsigned long [nMarker];
  Tag_to_Marker = new string [config->GetnMarker_Max()];

  /*--- Retrieve the name of the surface markers as well as
        the number of surface elements for every marker. ---*/
  const auto &sectionNames       = mesh->GetMarkerNames();
  const auto &nSurfElemPerMarker = mesh->GetNumberOfSurfaceElementsAllMarkers();

  /*--- Loop over all sections that we extracted from the grid file
        that were identified as boundary element sections so that we can
        store those elements into our SU2 data structures. ---*/
  for (int iMarker = 0; iMarker < nMarker; ++iMarker) {

    /*--- Get the string name and set the number of surface elements
          for this marker. ---*/
    string Marker_Tag    = sectionNames[iMarker];
    nElem_Bound[iMarker] = nSurfElemPerMarker[iMarker];

    /*--- Allocate the memory of the pointers for the surface
          elements for this marker. ---*/
    bound[iMarker] = new CPrimalGrid*[nElem_Bound[iMarker]];

    /*--- Retrieve the boundary element data for this marker. ---*/
    const auto &dataElems = mesh->GetSurfaceElementConnectivityForMarker(iMarker);

    /*--- Loop over the number of boundary elements for this marker. ---*/
    unsigned long ind = 0;
    for (unsigned long jElem=0; jElem<nElem_Bound[iMarker]; ++jElem) {

      /*--- Create a boundary FEM element from the data dataElems. ---*/
      const auto *dataElem = dataElems.data() + ind;
      bound[iMarker][jElem] = new CPrimalGridBoundFEM(dataElem);

      /*--- Update ind for the next element. ---*/
      ind += dataElem[2] + 5;
    }

    /*--- Update config file lists in order to store the boundary
          information for this marker in the correct place. ---*/
    Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag)] = Marker_Tag;
    config->SetMarker_All_TagBound(iMarker, Marker_Tag);
    config->SetMarker_All_KindBC(iMarker, config->GetMarker_CfgFile_KindBC(Marker_Tag));
    config->SetMarker_All_Monitoring(iMarker, config->GetMarker_CfgFile_Monitoring(Marker_Tag));
    config->SetMarker_All_GeoEval(iMarker, config->GetMarker_CfgFile_GeoEval(Marker_Tag));
    config->SetMarker_All_Designing(iMarker, config->GetMarker_CfgFile_Designing(Marker_Tag));
    config->SetMarker_All_Plotting(iMarker, config->GetMarker_CfgFile_Plotting(Marker_Tag));
    config->SetMarker_All_Analyze(iMarker, config->GetMarker_CfgFile_Analyze(Marker_Tag));
    config->SetMarker_All_ZoneInterface(iMarker, config->GetMarker_CfgFile_ZoneInterface(Marker_Tag));
    config->SetMarker_All_DV(iMarker, config->GetMarker_CfgFile_DV(Marker_Tag));
    config->SetMarker_All_Moving(iMarker, config->GetMarker_CfgFile_Moving(Marker_Tag));
    config->SetMarker_All_Deform_Mesh(iMarker, config->GetMarker_CfgFile_Deform_Mesh(Marker_Tag));
    config->SetMarker_All_Fluid_Load(iMarker, config->GetMarker_CfgFile_Fluid_Load(Marker_Tag));
    config->SetMarker_All_PyCustom(iMarker, config->GetMarker_CfgFile_PyCustom(Marker_Tag));
    config->SetMarker_All_PerBound(iMarker, config->GetMarker_CfgFile_PerBound(Marker_Tag));
    config->SetMarker_All_SendRecv(iMarker, NONE);
    config->SetMarker_All_Turbomachinery(iMarker, config->GetMarker_CfgFile_Turbomachinery(Marker_Tag));
    config->SetMarker_All_TurbomachineryFlag(iMarker, config->GetMarker_CfgFile_TurbomachineryFlag(Marker_Tag));
    config->SetMarker_All_MixingPlaneInterface(iMarker, config->GetMarker_CfgFile_MixingPlaneInterface(Marker_Tag));
  }
}

void CPhysicalGeometry::SetColorFEMGrid_Parallel(CConfig *config) {

  /*--- Initialize the color of the elements. ---*/
  for(unsigned long i=0; i<nElem; ++i)
    elem[i]->SetColor(0);

  /*--- Determine the standard elements for the grid. ---*/
  DetermineFEMStandardElements(config);

  /*--- Determine the matching faces of the local elements. ---*/
  vector<CFaceOfElement> localMatchingFaces;
  DetermineMatchingFacesFEMGrid(config, localMatchingFaces);

  /*--- In the function above the periodic boundaries are not found.
        A different treatment must be used in order to find these. ---*/
  DeterminePeriodicFacesFEMGrid(config, localMatchingFaces);

  /*--- So far only the matching faces have been found. Now find the
        non-matching faces of the local elements. ---*/
  DetermineNonMatchingFacesFEMGrid(config, localMatchingFaces);

  /*--- Determine whether or not the Jacobians of the elements and faces
        can be considered constant. ---*/
  DetermineFEMConstantJacobiansAndLenScale(config);

  /*--- Determine the donor elements for the wall function treatment. ---*/
  DetermineDonorElementsWallFunctions(config);

  /*--- Determine the time levels of the elements. This is only relevant
        when time accurate local time stepping is employed. ---*/
  map<unsigned long, CUnsignedShort2T> mapExternalElemIDToTimeLevel;
  DetermineTimeLevelElements(config, localMatchingFaces, mapExternalElemIDToTimeLevel);

  /*--- Determine the ownership of the internal faces, i.e. which adjacent
        element is responsible for computing the fluxes through the face. ---*/
  DetermineOwnershipInternalFaces(localMatchingFaces, mapExternalElemIDToTimeLevel);

  /*--- All the matching face information is known now, including periodic
        faces. Store the information of the neighbors in the data structure
        for the local elements. ---*/
  StoreFaceInfoInLocalElements(localMatchingFaces);

  /*--- Create the vector of vectors that describe the connectivity
        of the graph. ---*/
  vector<vector<unsigned long> > adjacency;
  DetermineGraphAdjacency(localMatchingFaces, adjacency);

  /*--- Determine the weigts of the graph. ---*/
  vector<passivedouble> vwgt;
  vector<vector<passivedouble> > adjwgt;
  DetermineFEMGraphWeights(config, localMatchingFaces, adjacency,
                           mapExternalElemIDToTimeLevel, vwgt, adjwgt);

  /*--- The remainder of this function should only be called if we have
        parallel support with MPI. ---*/
#ifdef HAVE_MPI

  /*--- If the ParMETIS library is compiled and linked, the colors
        or determined via ParMETIS. Otherwise a linear distribution
        is used, which is usually inefficient. ---*/
#ifdef HAVE_PARMETIS
  DetermineFEMColorsViaParMETIS(adjacency, vwgt, adjwgt);
#else
  if(size > SINGLE_NODE)
  {
    if (rank == MASTER_NODE) {
      cout << endl;
      cout << "--------------------- WARNING -------------------------------" << endl;
      cout << "SU2 compiled without PARMETIS. A linear distribution is used." << endl;
      cout << "This is very inefficient" << endl;
      cout << "-------------------------------------------------------------" << endl;
      cout << endl;
    }

    for(unsigned long i=0; i<nElem; ++i)
      elem[i]->SetColor(rank);
  }

#endif  /* HAVE_PARMETIS */
#endif  /* HAVE_MPI */
}

void CPhysicalGeometry::AllocateMemoryMatricesMetrics(
                                vector<CFEMStandardElementBase *>           &standardElements,
                                const unsigned short                        nDer,
                                vector<su2activevector>                     *matricesJacobians,
                                vector<ColMajorMatrix<su2double> >          *matricesCoor,
                                vector<ColMajorMatrix<su2double> >          *matricesNormals,
                                vector<vector<ColMajorMatrix<su2double> > > *matricesDerCoor) {

  /*--- Allocate the memory for matricesJacobians, if needed. ---*/
  if( matricesJacobians ) {

    matricesJacobians->resize(standardElements.size());

    for(unsigned long i=0; i<standardElements.size(); ++i) {
      const unsigned short sizeInt = standardElements[i]->GetNIntegrationPad();
      matricesJacobians->at(i).resize(sizeInt);
      matricesJacobians->at(i).setConstant(0.0);
    }
  }

  /*--- Allocate the memory for matricesCoor, if needed. ---*/
  if( matricesCoor ) {

    matricesCoor->resize(standardElements.size());

    for(unsigned long i=0; i<standardElements.size(); ++i) {
      const unsigned short sizeDOF = standardElements[i]->GetNDOFs();
      matricesCoor->at(i).resize(sizeDOF, nDim);
      matricesCoor->at(i).setConstant(0.0);
    }
  }

  /*--- Allocate the memory for matricesNormals, if needed. ---*/
  if( matricesNormals ) {

    matricesNormals->resize(standardElements.size());

    for(unsigned long i=0; i<standardElements.size(); ++i) {
      const unsigned short sizeInt = standardElements[i]->GetNIntegrationPad();
      matricesNormals->at(i).resize(sizeInt, nDim);
      matricesNormals->at(i).setConstant(0.0);
    }
  }

  /*--- Allocate the memory for matricesDerCoor, if needed. ---*/
  if( matricesDerCoor ) {

    matricesDerCoor->resize(standardElements.size());

    for(unsigned long i=0; i<standardElements.size(); ++i) {
      const unsigned short sizeInt = standardElements[i]->GetNIntegrationPad();
      matricesDerCoor->at(i).resize(nDer);

      for(unsigned short j=0; j<nDer; ++j) { 
        matricesDerCoor->at(i)[j].resize(sizeInt, nDim);

        /*--- Set the default values to avoid problems later on. ---*/
        matricesDerCoor->at(i)[j].setConstant(0.0);

        for(unsigned short k=0; k<sizeInt; ++k)
          matricesDerCoor->at(i)[j](k,j) = 1.0;
      }
    }
  }
}

void CPhysicalGeometry::DetermineDonorElementsWallFunctions(CConfig *config) {

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Check whether wall functions are used at all.              ---*/
  /*--------------------------------------------------------------------------*/

  unsigned long nWallFaces = 0;
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

    switch (config->GetMarker_All_KindBC(iMarker)) {
      case ISOTHERMAL:
      case HEAT_FLUX: {
        const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if(config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::NONE)
          nWallFaces += nElem_Bound[iMarker];
        break;
      }
      default:  /*--- Just to avoid a compiler warning. ---*/
        break;
    }
  }

  /*--- If no wall functions are used, nothing needs to be done and a
        return can be made. ---*/
  if(nWallFaces == 0) return;

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Build the ADT of the linear sub-elements of the locally    ---*/
  /*---         stored volume elements.                                    ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Determine a mapping from the global point ID to the local index
        of the points. ---*/
  map<unsigned long,unsigned long> globalPointIDToLocalInd;
  for(unsigned long i=0; i<nPoint; ++i) {
    globalPointIDToLocalInd[nodes->GetGlobalIndex(i)] = i;
  }

  /*--- Define the vectors, which store the mapping from the subelement to the
        parent element, subelement ID within the parent element, the element
        type and the connectivity of the subelements. ---*/
  vector<unsigned long>  parentElement;
  vector<unsigned short> subElementIDInParent;
  vector<unsigned short> VTK_TypeElem;
  vector<unsigned long>  elemConn;

  /*--- Loop over the local volume elements to create the connectivity of
        the linear sub-elements. ---*/
  for(unsigned long l=0; l<nElem; ++l) {

    /*--- Determine the standard element of the volume element. ---*/
    const unsigned short VTK_Parent = elem[l]->GetVTK_Type();
    const unsigned short nPolyGrid  = elem[l]->GetNPolyGrid();
    const unsigned short orderExact = config->GetOrderExactIntegrationFEM(nPolyGrid, false);

    unsigned long ii;
    for(ii=0; ii<standardVolumeElements.size(); ++ii)
      if( standardVolumeElements[ii]->SameStandardElement(VTK_Parent, nPolyGrid, orderExact) )
        break;

    /*-- Determine the necessary data for splitting the element in its linear
         sub-elements. ---*/
    unsigned short VTK_Type[]  = {standardVolumeElements[ii]->GetVTK_SubType1(),
                                  standardVolumeElements[ii]->GetVTK_SubType2()};
    unsigned short nSubElems[] = {0, 0};
    unsigned short nDOFsPerSubElem[] = {0, 0};
    const unsigned short *connSubElems[] = {nullptr, nullptr};

    if(VTK_Type[0] != NONE) {
      nSubElems[0]       = standardVolumeElements[ii]->GetNSubElemsType1();
      nDOFsPerSubElem[0] = standardVolumeElements[ii]->GetNDOFsPerSubElem(VTK_Type[0]);
      connSubElems[0]    = standardVolumeElements[ii]->GetSubConnType1();
    }

    if(VTK_Type[1] != NONE) {
      nSubElems[1]       = standardVolumeElements[ii]->GetNSubElemsType2();
      nDOFsPerSubElem[1] = standardVolumeElements[ii]->GetNDOFsPerSubElem(VTK_Type[1]);
      connSubElems[1]    = standardVolumeElements[ii]->GetSubConnType2();
    }

    /*--- Store the connectivity of the sub-elements. Note that local node
          numbering must be used for these sub-elements. ---*/
    unsigned short jj = 0;
    for(unsigned short i=0; i<2; ++i) {
      unsigned short kk = 0;
      for(unsigned short j=0; j<nSubElems[i]; ++j, ++jj) {
        parentElement.push_back(elem[l]->GetGlobalElemID());
        subElementIDInParent.push_back(jj);
        VTK_TypeElem.push_back(VTK_Type[i]);

        for(unsigned short k=0; k<nDOFsPerSubElem[i]; ++k, ++kk) {
          unsigned long nodeID = elem[l]->GetNode(connSubElems[i][kk]);
          map<unsigned long,unsigned long>::const_iterator MI;
          MI = globalPointIDToLocalInd.find(nodeID);

          elemConn.push_back(MI->second);
        }
      }
    }
  }

  /*--- Store the coordinates of the locally stored nodes in the format
        expected by the ADT. ---*/
  vector<su2double> volCoor(nDim*nPoint);

  unsigned long jj = 0;
  for(unsigned long l=0; l<nPoint; ++l) {
    for(unsigned short k=0; k<nDim; ++k, ++jj)
      volCoor[jj] = nodes->GetCoord(l, k);
  }

  /*--- Build the local ADT. ---*/
  CADTElemClass localVolumeADT(nDim, volCoor, elemConn, VTK_TypeElem,
                               subElementIDInParent, parentElement, false);

  /*--- Release the memory of the vectors used to build the ADT. To make sure
        that all the memory is deleted, the swap function is used. ---*/
  vector<unsigned short>().swap(subElementIDInParent);
  vector<unsigned short>().swap(VTK_TypeElem);
  vector<unsigned long>().swap(parentElement);
  vector<unsigned long>().swap(elemConn);
  vector<su2double>().swap(volCoor);

  /*--------------------------------------------------------------------------*/
  /*--- Step 3. Search for donor elements at the exchange locations in     ---*/
  /*---         the local elements.                                        ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Determine whether or not the LGL node distribution is used. ---*/
  const bool useLGL = config->GetKind_FEM_GridDOFsLocation() == LGL;

  /*--- Define the vectors, which store the boundary marker, boundary element ID
        and exchange coordinates for the integration points for which no donor
        element was found in the locally stored volume elements. ---*/
  vector<unsigned short> markerIDGlobalSearch;
  vector<unsigned long>  boundaryElemIDGlobalSearch;
  vector<su2double>      coorExGlobalSearch;

  /*---- Start of the OpenMP parallel region, if supported. ---*/
  SU2_OMP_PARALLEL
  {
    /*--- Define the matrices used to store the coordinates, its derivatives,
          the Jacobians and the unit normals. Every standard element gets its
          own set of matrices, such that the performance is maximized. Not so
          important here, but it is compatible with the approach used in the
          computationally intensive part. ---*/
    vector<su2activevector> matricesJacobians;
    vector<ColMajorMatrix<su2double> > matricesCoor, matricesNormals;
    vector<vector<ColMajorMatrix<su2double> > > matricesDerCoor;

    /*--- Allocate the memory for the matrices to compute the face metrics. ---*/
    AllocateMemoryMatricesMetrics(standardFaceElements, nDim-1, &matricesJacobians,
                                  &matricesCoor, &matricesNormals, &matricesDerCoor);

    /*--- Loop over the markers and select the ones for which a wall function
          treatment must be carried out. ---*/
    for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

      switch (config->GetMarker_All_KindBC(iMarker)) {
        case ISOTHERMAL:
        case HEAT_FLUX: {
          const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if(config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::NONE) {

            /*--- Retrieve the floating point information for this boundary marker.
                  The exchange location is the first element of this array. ---*/
            const su2double *doubleInfo = config->GetWallFunction_DoubleInfo(Marker_Tag);

            /*--- Determine the chunk size for the OMP loop below. ---*/
#ifdef HAVE_OMP
            const size_t omp_chunk_size = computeStaticChunkSize(nElem_Bound[iMarker],
                                                                 omp_get_num_threads(), 64);
#endif
            /*--- Loop over the local boundary elements for this marker. ---*/
            SU2_OMP_FOR_DYN(omp_chunk_size)
            for(unsigned long l=0; l<nElem_Bound[iMarker]; ++l) {

              /*--- Determine the local volume element of this boundary element. ---*/
              const unsigned long globalElemID = bound[iMarker][l]->GetDomainElement();
              unordered_map<unsigned long,unsigned long>::const_iterator MI;
              MI = Global_to_Local_Elem.find(globalElemID);

              if(MI == Global_to_Local_Elem.end())
                SU2_MPI::Error(string("Global element ID not found. This should not happen"),
                               CURRENT_FUNCTION);

              const unsigned long elemID = MI->second;

              /*--- Get the corner points of the boundary element. Note that this
                    is an overloaded function, hence the arguments allow for
                    multiple faces ---. */
              unsigned short nFaces;
              unsigned short nPointsPerFace[6];
              unsigned long  faceConn[6][4];
              bound[iMarker][l]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

              /*--- Create an object of CFaceOfElement to store the information. ---*/
              CFaceOfElement boundaryFace;
              boundaryFace.nCornerPoints = nPointsPerFace[0];
              for(unsigned short i=0; i<nPointsPerFace[0]; ++i)
                boundaryFace.cornerPoints[i] = faceConn[0][i];
              boundaryFace.elemID0 = elemID;

              /*--- Renumber the corner points of the face, but keep the orientation. ---*/
              boundaryFace.CreateUniqueNumberingWithOrientation();
              const bool boundaryFaceSwapped = boundaryFace.elemID1 == elemID;

              /*--- Get the corner points of all the faces from the corresponding
                    volume element. ---*/
              elem[elemID]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

              /*--- Loop over the faces of the element to determine which one is
                    the boundary face. ---*/
              unsigned long ii;
              bool outwardPointing = true;   /* Initialized to avoid compiler warning. */
              for(ii=0; ii<nFaces; ++ii) {

                /*--- Create an object of CFaceOfElement to store the information.
                      Renumber the corner points, but keep the orientation. ---*/
                CFaceOfElement thisFace;
                thisFace.nCornerPoints = nPointsPerFace[ii];
                for(unsigned short i=0; i<nPointsPerFace[ii]; ++i)
                  thisFace.cornerPoints[i] = faceConn[ii][i];
                thisFace.elemID0 = elemID;

                thisFace.CreateUniqueNumberingWithOrientation();
                const bool thisFaceSwapped = thisFace.elemID1 == elemID;

                /*--- Check if this face is the boundary face. ---*/
                if(boundaryFace == thisFace) {

                  /*--- Determine whether the orientation of the boundary face is
                        pointing out of the adjacent element. ---*/
                  outwardPointing = boundaryFaceSwapped == thisFaceSwapped;

                  /* Break the loop over the faces of the element. */
                  break;
                }
              }

              /*--- Additional check, just to be sure. ---*/
              if(ii == nFaces)
                SU2_MPI::Error(string("Boundary face not found in faces of element. This should not happen"),
                               CURRENT_FUNCTION);

              /*--- Determine whether or not the Jacobian of the boundary face
                    is constant. ---*/
              const bool constJac = elem[elemID]->GetJacobianConstantFace(ii);
              bound[iMarker][l]->SetJacobianConsideredConstant(constJac);

              /*--- Abbreviate the boundary element type, polynomial degree of the
                    solution and determine the polynomial degree that must be integrated
                    exactly, based on the polynomial degree of the solution. ---*/
              const unsigned short VTK_Type   = bound[iMarker][l]->GetVTK_Type();
              const unsigned short nPolySol   = elem[elemID]->GetNPolySol();
              const unsigned short orderExact = config->GetOrderExactIntegrationFEM(nPolySol, constJac);

              /*--- Determine the corresponding standard element of the face. ---*/
              unsigned long jj;
              for(jj=0; jj<standardFaceElements.size(); ++jj)
                if( standardFaceElements[jj]->SameStandardElement(VTK_Type, nPolySol, orderExact) )
                  break;

              /*--- Store the coordinates of the face in the correct entry of matricesCoor. ---*/
              const unsigned short nDOFs = standardFaceElements[jj]->GetNDOFs();
              for(unsigned short j=0; j<nDOFs; ++j) {
                unsigned long nodeID = bound[iMarker][l]->GetNode(j);
                map<unsigned long,unsigned long>::const_iterator MI;
                MI = globalPointIDToLocalInd.find(nodeID);
                nodeID = MI->second;
                for(unsigned short k=0; k<nDim; ++k)
                  matricesCoor[jj](j,k) = nodes->GetCoord(nodeID, k);
              }

              /*--- Create the coordinates and the unit normals in the integration points.
                    Make sure to compute the normals first, because then the memory of
                    matricesDerCoor can be used to store the coordinates of the
                    integration points. Also store the unit normals a bit easier. ---*/
              standardFaceElements[jj]->UnitFaceNormals(useLGL, matricesCoor[jj],
                                                        matricesDerCoor[jj], matricesNormals[jj],
                                                        matricesJacobians[jj]);

              ColMajorMatrix<su2double> &coorInt  = matricesDerCoor[jj][0];
              ColMajorMatrix<su2double> &unitNorm = matricesNormals[jj];
              standardFaceElements[jj]->CoorIntPoints(useLGL, matricesCoor[jj], coorInt);

              /*--- Set the multiplication factor for the normal, such that
                    an inward pointing normal is obtained. It is scaled with
                    the exchange distance. ---*/
              const su2double factNorm = outwardPointing ? -doubleInfo[0] : doubleInfo[0];

              /*--- Determine the number of integration points and reserve the
                    memory to store the donor elements of the face. ---*/
              const unsigned short nInt = standardFaceElements[jj]->GetNIntegration();
              vector<unsigned long> donorElementsFace;
              donorElementsFace.reserve(nInt);

              /*--- Loop over the number of integration points. ---*/
              for(unsigned short i=0; i<nInt; ++i) {

                /*--- Determine the coordinate of the exchange point. ---*/
                su2double coorEx[3] = {0.0, 0.0, 0.0};
                for(unsigned short k=0; k<nDim; ++k)
                  coorEx[k] = coorInt(i,k) + factNorm*unitNorm(i,k);

                /*--- Search for the element, which contains the exchange location. ---*/
                unsigned short subElem;
                unsigned long  parElem;
                int            mpirank;
                su2double      parCoor[3], weightsInterpol[8];

                if( localVolumeADT.DetermineContainingElement(coorEx, subElem,
                                                              parElem, mpirank, parCoor,
                                                              weightsInterpol) ) {

                  /*--- Donor element found. Store it in donorElementsFace. ---*/
                  donorElementsFace.push_back(parElem);
                }
                else {

                  /*--- Donor element not found in the local ADT. Store the exchange
                        coordinates, boundary marker and boundary element ID. ---*/
                  SU2_OMP_CRITICAL
                  {
                    markerIDGlobalSearch.push_back(iMarker);
                    boundaryElemIDGlobalSearch.push_back(l);
                    for(unsigned short iDim=0; iDim<nDim; ++iDim)
                      coorExGlobalSearch.push_back(coorEx[iDim]);
                  }
		  END_SU2_OMP_CRITICAL
                }
              }

              /*--- Sort donorElementsFace in increasing order and remove the
                    the double entities. ---*/
              sort(donorElementsFace.begin(), donorElementsFace.end());
              vector<unsigned long>::iterator lastEntry;
              lastEntry = unique(donorElementsFace.begin(), donorElementsFace.end());
              donorElementsFace.erase(lastEntry, donorElementsFace.end()); 

              /*--- Store the donor elements in the data structure for
                    this boundary element. ---*/
              bound[iMarker][l]->SetDonorsWallFunctions(donorElementsFace);
            }
	    END_SU2_OMP_FOR
          }

          break;
        }

        default:  /*--- Just to avoid a compiler warning. ---*/
          break;
      }
    }

  }
  END_SU2_OMP_PARALLEL

  /*--- The remaining part of this function only needs to be carried out in parallel mode. ---*/
#ifdef HAVE_MPI

  /*--- Determine the number of search points for which a global search must be
        carried out for each rank and store them in such a way that the info can
        be used directly in Allgatherv. ---*/
  vector<int> recvCounts(size), displs(size);
  int nLocalSearchPoints = static_cast<int>(markerIDGlobalSearch.size());

  SU2_MPI::Allgather(&nLocalSearchPoints, 1, MPI_INT, recvCounts.data(), 1,
                     MPI_INT, SU2_MPI::GetComm());
  displs[0] = 0;
  for(int i=1; i<size; ++i) displs[i] = displs[i-1] + recvCounts[i-1];

  const int nGlobalSearchPoints = displs.back() + recvCounts.back();

  /*--- Check if there actually are global searches to be carried out. ---*/
  if(nGlobalSearchPoints > 0) {

    /*--- Create a cumulative storage version of recvCounts. ---*/
    vector<int> nSearchPerRank(size+1);
    nSearchPerRank[0] = 0;

    for(int i=0; i<size; ++i)
      nSearchPerRank[i+1] = nSearchPerRank[i] + recvCounts[i];

    /*--- Gather the data of the search points for which a global search must
          be carried out on all ranks. ---*/
    vector<unsigned short> bufMarkerIDGlobalSearch(nGlobalSearchPoints);
    SU2_MPI::Allgatherv(markerIDGlobalSearch.data(), nLocalSearchPoints,
                        MPI_UNSIGNED_SHORT, bufMarkerIDGlobalSearch.data(),
                        recvCounts.data(), displs.data(), MPI_UNSIGNED_SHORT,
                        SU2_MPI::GetComm());

    vector<unsigned long> bufBoundaryElemIDGlobalSearch(nGlobalSearchPoints);
    SU2_MPI::Allgatherv(boundaryElemIDGlobalSearch.data(), nLocalSearchPoints,
                        MPI_UNSIGNED_LONG, bufBoundaryElemIDGlobalSearch.data(),
                        recvCounts.data(), displs.data(), MPI_UNSIGNED_LONG,
                        SU2_MPI::GetComm());

    for(int i=0; i<size; ++i) {recvCounts[i] *= nDim; displs[i] *= nDim;}
    vector<su2double> bufCoorExGlobalSearch(nDim*nGlobalSearchPoints);
    SU2_MPI::Allgatherv(coorExGlobalSearch.data(), nDim*nLocalSearchPoints,
                        MPI_DOUBLE, bufCoorExGlobalSearch.data(),
                        recvCounts.data(), displs.data(), MPI_DOUBLE,
                        SU2_MPI::GetComm());

    /*--- Buffers to store the return information. ---*/
    vector<unsigned short> markerIDReturn;
    vector<unsigned long>  boundaryElemIDReturn;
    vector<unsigned long>  volElemIDDonorReturn;

    /*--- Loop over the number of global search points to check if these points
          are contained in the volume elements of this rank. The loop is carried
          out as a double loop, such that the rank where the point resides is
          known as well. Furthermore, it is not necessary to search the points
          that were not found earlier on this rank. The vector recvCounts is used
          as storage for the number of search items that must be returned to the
          other ranks. ---*/
    for(int rankID=0; rankID<size; ++rankID) {
      recvCounts[rankID] = 0;
      if(rankID != rank) {
        for(int i=nSearchPerRank[rankID]; i<nSearchPerRank[rankID+1]; ++i) {

          /*--- Search the local ADT for the coordinate of the exchange point
                and check if it is found. ---*/
          unsigned short subElem;
          unsigned long  parElem;
          int            rankDonor;
          su2double      parCoor[3], weightsInterpol[8];
          if( localVolumeADT.DetermineContainingElement(bufCoorExGlobalSearch.data() + i*nDim,
                                                        subElem, parElem, rankDonor, parCoor,
                                                        weightsInterpol) ) {

            /*--- Store the required data in the return buffers. ---*/
            ++recvCounts[rankID];
            markerIDReturn.push_back(bufMarkerIDGlobalSearch[i]);
            boundaryElemIDReturn.push_back(bufBoundaryElemIDGlobalSearch[i]);
            volElemIDDonorReturn.push_back(parElem);
          }
        }
      }
    }

    /*--- Create a cumulative version of recvCounts. ---*/
    for(int i=0; i<size; ++i)
      nSearchPerRank[i+1] = nSearchPerRank[i] + recvCounts[i];

    /*--- Determine the number of return messages this rank has to receive.
          Use displs and recvCounts as temporary storage. ---*/
    int nRankSend = 0;
    for(int i=0; i<size; ++i) {
      if( recvCounts[i] ) {recvCounts[i] = 1; ++nRankSend;}
      displs[i] = 1;
    }

    int nRankRecv;
    SU2_MPI::Reduce_scatter(recvCounts.data(), &nRankRecv, displs.data(),
                            MPI_INT, MPI_SUM, SU2_MPI::GetComm());

    /*--- Send the data using nonblocking sends to avoid deadlock. ---*/
    vector<SU2_MPI::Request> commReqs(3*nRankSend);
    nRankSend = 0;
    for(int i=0; i<size; ++i) {
      if( recvCounts[i] ) {
        const int sizeMessage = nSearchPerRank[i+1] - nSearchPerRank[i];
        SU2_MPI::Isend(markerIDReturn.data() + nSearchPerRank[i],
                       sizeMessage, MPI_UNSIGNED_SHORT, i, i, SU2_MPI::GetComm(),
                       &commReqs[nRankSend++]);
        SU2_MPI::Isend(boundaryElemIDReturn.data() + nSearchPerRank[i],
                       sizeMessage, MPI_UNSIGNED_LONG, i, i+1, SU2_MPI::GetComm(),
                       &commReqs[nRankSend++]);
        SU2_MPI::Isend(volElemIDDonorReturn.data() + nSearchPerRank[i],
                       sizeMessage, MPI_UNSIGNED_LONG, i, i+2, SU2_MPI::GetComm(),
                       &commReqs[nRankSend++]);
      }
    }

    /*--- Loop over the number of ranks from which I receive return data. ---*/
    for(int i=0; i<nRankRecv; ++i) {

      /*--- Block until a message with unsigned shorts arrives from any processor.
            Determine the source and the size of the message. ---*/
      SU2_MPI::Status status;
      SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
      int source = status.MPI_SOURCE;

      int sizeMess;
      SU2_MPI::Get_count(&status, MPI_UNSIGNED_SHORT, &sizeMess);

      /*--- Allocate the memory for the receive buffers. ---*/
      vector<unsigned short> bufMarkerIDReturn(sizeMess);
      vector<unsigned long>  bufBoundaryElemIDReturn(sizeMess);
      vector<unsigned long>  bufVolElemIDDonorReturn(sizeMess);

      /*--- Receive the three messages using blocking receives. ---*/
      SU2_MPI::Recv(bufMarkerIDReturn.data(), sizeMess, MPI_UNSIGNED_SHORT,
                    source, rank, SU2_MPI::GetComm(), &status);

      SU2_MPI::Recv(bufBoundaryElemIDReturn.data(), sizeMess, MPI_UNSIGNED_LONG,
                    source, rank+1, SU2_MPI::GetComm(), &status);

      SU2_MPI::Recv(bufVolElemIDDonorReturn.data(), sizeMess, MPI_UNSIGNED_LONG,
                    source, rank+2, SU2_MPI::GetComm(), &status);

      /*--- Loop over the data just received and add it to the wall function
            donor information of the corresponding boundary element. ---*/
      for(int j=0; j<sizeMess; ++j) {
        const unsigned short iMarker = bufMarkerIDReturn[j];
        const unsigned long  l       = bufBoundaryElemIDReturn[j];
        const unsigned long  volID   = bufVolElemIDDonorReturn[j];

        bound[iMarker][l]->AddDonorWallFunctions(volID);
      }
    }

    /*--- Complete the non-blocking sends. ---*/
    SU2_MPI::Waitall(nRankSend, commReqs.data(), MPI_STATUSES_IGNORE);

    /*--- Wild cards have been used in the communication,
          so synchronize the ranks to avoid problems. ---*/
    SU2_MPI::Barrier(SU2_MPI::GetComm());

    /*--- Loop again over the boundary elements of the marker for which a wall
          function treatment must be used and make remove the multiple entries
          of the donor information. ---*/
    for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

      switch (config->GetMarker_All_KindBC(iMarker)) {
        case ISOTHERMAL:
        case HEAT_FLUX: {
          const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if(config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::NONE) {

            for(unsigned long l=0; l<nElem_Bound[iMarker]; ++l)
              bound[iMarker][l]->RemoveMultipleDonorsWallFunctions();
          }

          break;
        }

        default:  /*--- Just to avoid a compiler warning. ---*/
          break;
      }
    }
  }

#endif   /*--- HAVE_MPI ---*/

}

#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS
void CPhysicalGeometry::DetermineFEMColorsViaParMETIS(vector<vector<unsigned long> > &adjacency,
                                                      vector<passivedouble>          &vwgt,
                                                      vector<vector<passivedouble> > &adjwgt) {

  /*--- Only call ParMETIS if we have more than one rank to avoid errors ---*/
  if(size > SINGLE_NODE) {

    /*--- Define the linear partitioning of the elements. ---*/
    CLinearPartitioner elemPartitioner(Global_nElem, 0);

    /*--- Determine the array, which stores the distribution of the graph nodes
          over the ranks.     ---*/
    vector<idx_t> vtxdist(size+1);
    for(int i=0; i<=size; ++i)
      vtxdist[i] = static_cast<idx_t>(elemPartitioner.GetCumulativeSizeBeforeRank(i));

    /*--- Create the array xadjPar, which contains the number of edges for each
          vertex of the graph in ParMETIS format. ---*/
    vector<idx_t> xadjPar(nElem+1);
    xadjPar[0] = 0;
    for(unsigned long i=0; i<nElem; ++i)
      xadjPar[i+1] = xadjPar[i] + static_cast<idx_t>(adjacency[i].size());

    /*--- Create the adjacency information in ParMETIS format. ---*/
    vector<idx_t> adjacencyPar(xadjPar[nElem]);
    unsigned long ii = 0;
    for(unsigned long i=0; i<nElem; ++i) {
      for(unsigned long j=0; j<adjacency[i].size(); ++j, ++ii)
        adjacencyPar[ii] = static_cast<idx_t>(adjacency[i][j]);
    }

    /*--- Create the vertex weights in ParMETIS format. ---*/
    vector<idx_t> vwgtPar(vwgt.size());
    for(unsigned long i=0; i<vwgt.size(); ++i)
      vwgtPar[i] = static_cast<idx_t>(ceil(vwgt[i]));

    /*--- Create the adjacency weight in ParMETIS format. ---*/
    vector<idx_t> adjwgtPar(xadjPar[nElem]);
    ii = 0;
    for(unsigned long i=0; i<nElem; ++i) {
      for(unsigned long j=0; j<adjwgt[i].size(); ++j, ++ii)
        adjwgtPar[ii] = static_cast<idx_t>(ceil(adjwgt[i][j]));
    }

    /*--- The scalar variables and the options array for the call to ParMETIS. ---*/
    idx_t  wgtflag = 3;                         // Weights on both the vertices and edges.
    idx_t  numflag = 0;                         // C-numbering.
    idx_t  ncon    = 2;                         // Number of constraints.
    real_t ubvec[] = {1.05, 1.05};              // Tolerances for the vertex weights, recommended value is 1.05.
    idx_t  nparts  = static_cast<idx_t>(size);  // Number of subdomains. Must be number of MPI ranks.
    idx_t  options[METIS_NOPTIONS];             // Just use the default options.
    METIS_SetDefaultOptions(options);
    options[1] = 0;

    /*--- Make sure that an equal distribution is obtained. ---*/
    vector<real_t> tpwgts(size*ncon, 1.0/(static_cast<real_t>(size)));

    /*--- Calling ParMETIS ---*/
    vector<idx_t> part(nElem);
    if (rank == MASTER_NODE) cout << "Calling ParMETIS...";

    idx_t edgecut;
    MPI_Comm comm = SU2_MPI::GetComm();
    ParMETIS_V3_PartKway(vtxdist.data(), xadjPar.data(), adjacencyPar.data(),
                         vwgtPar.data(), adjwgtPar.data(), &wgtflag, &numflag,
                         &ncon, &nparts, tpwgts.data(), ubvec, options,
                         &edgecut, part.data(), &comm);
    if (rank == MASTER_NODE) {
      cout << " graph partitioning complete (";
      cout << edgecut << " edge cuts)." << endl;
    }

    /*--- Set the color of the elements to the outcome of ParMETIS. ---*/
    for(unsigned long i=0; i<nElem; ++i)
      elem[i]->SetColor(part[i]);
  }
}
#endif
#endif

void CPhysicalGeometry::DetermineFEMConstantJacobiansAndLenScale(CConfig *config) {

  /*--- Determine a mapping from the global point ID to the local index
        of the points.    ---*/
  map<unsigned long,unsigned long> globalPointIDToLocalInd;
  for(unsigned long i=0; i<nPoint; ++i) {
    globalPointIDToLocalInd[nodes->GetGlobalIndex(i)] = i;
  }

  /*--- Determine the chunk size for the OMP loop below. ---*/
#ifdef HAVE_OMP
  const size_t omp_chunk_size = computeStaticChunkSize(nElem, omp_get_num_threads(), 64);
#endif

  /*--- Vectors to store the minimum and maximum Jacobian of the elements,
        both for the LGL and equidistant node distribution. ---*/
  vector<su2double> jacMinLGL(nElem),  jacMaxLGL(nElem);
  vector<su2double> jacMinEqui(nElem), jacMaxEqui(nElem);

  /*--- Counters, which keep track of the number of different grid
        location types present. Initialized to zero. ---*/
  unsigned long counterGridLocation[NO_PREFERRED_LOCATION+1] = {0};

  /*---- Start of the OpenMP parallel region, if supported. ---*/
  SU2_OMP_PARALLEL
  {
    /*--- Definition of the thread local version of counterGridLocation to handle
          the reduction handled manually to avoid complications with CODIPACK. ---*/
    unsigned long counterLocation[NO_PREFERRED_LOCATION+1] = {0};

    /*--- Define the matrices used to store the coordinates, its derivatives and
          the Jacobians. For later purposes, when computing the face metrics,
          als the matrices for storing the normals is defined. Every standard
          element get its own set of matrices, such that the performance is
          maximized. Not so important here, but it is compatible with the
          approach used in the computationally intensive part. ---*/
    vector<su2activevector> matricesJacobians;
    vector<ColMajorMatrix<su2double> > matricesCoor, matricesNormals;
    vector<vector<ColMajorMatrix<su2double> > > matricesDerCoor;

    /*--- Allocate the memory for these matrices. ---*/
    AllocateMemoryMatricesMetrics(standardVolumeElements, nDim, &matricesJacobians,
                                  &matricesCoor, nullptr, &matricesDerCoor);

    /*--- Loop over the local volume elements to compute the minimum and maximum
          Jacobian for both the LGL and equidistant point distribution and to
          determine the number of different location types. ---*/
    SU2_OMP_FOR_DYN(omp_chunk_size)
    for(unsigned long i=0; i<nElem; ++i) {

      /*--- Determine the order of the polynomials that must be integrated exactly.
            As the intention of this function is to find out whether or not the
            Jacobians are constant, it suffices to look at the polynomial degree
            of the grid. ---*/
      const unsigned short nPolyGrid  = elem[i]->GetNPolyGrid();
      const unsigned short orderExact = config->GetOrderExactIntegrationFEM(nPolyGrid, false);

      /*--- Determine the standard element of the volume element. ---*/
      unsigned long ii;
      for(ii=0; ii<standardVolumeElements.size(); ++ii)
        if( standardVolumeElements[ii]->SameStandardElement(elem[i]->GetVTK_Type(),
                                                            nPolyGrid, orderExact) )
          break;

      /*--- Retrieve the number of grid DOFs for this element. ---*/
      const unsigned short nDOFs = standardVolumeElements[ii]->GetNDOFs();

      /*--- Copy the coordinates into matricesCoor[ii], such that the gemm
            routines can be used to compute the derivatives. ---*/
      for(unsigned short j=0; j<nDOFs; ++j) {
        const unsigned long nodeID = elem[i]->GetNode(j);

        map<unsigned long,unsigned long>::const_iterator MI = globalPointIDToLocalInd.find(nodeID);
        const unsigned long ind = MI->second;
        for(unsigned short k=0; k<nDim; ++k)
          matricesCoor[ii](j,k) = nodes->GetCoord(ind, k);
      }

      /*--- Determine the minimum and maximum values of the Jacobians of the
            transformation to the standard element for the LGL and equidistant
            distribution of the grid DOFs. ---*/
      standardVolumeElements[ii]->MinMaxJacobians(true, matricesCoor[ii],
                                                  matricesDerCoor[ii],
                                                  matricesJacobians[ii],
                                                  jacMinLGL[i], jacMaxLGL[i]);

      standardVolumeElements[ii]->MinMaxJacobians(false, matricesCoor[ii],
                                                  matricesDerCoor[ii],
                                                  matricesJacobians[ii],
                                                  jacMinEqui[i], jacMaxEqui[i]);

      /*--- Determine the situation for the grid location and update
            the appropriate entry in counterLocation. ---*/
      if((jacMinLGL[i] <= 0.0) && (jacMinEqui[i] <= 0.0)) ++counterLocation[NO_VALID_LOCATION];
      else if(jacMinEqui[i] <= 0.0)                       ++counterLocation[LGL_ONLY];
      else if(jacMinLGL[i]  <= 0.0)                       ++counterLocation[EQUI_ONLY];
      else {

        /*--- Both point distribution produce a valid mapping. Pick the
              preferred one based on the ratio of the minimum and
              maximum Jacobian that occurs. ---*/
        const su2double ratioJacLGL = jacMaxLGL[i] /jacMinLGL[i];
        const su2double ratioJacEq  = jacMaxEqui[i]/jacMinEqui[i];

        if(     ratioJacEq/ratioJacLGL >= 1.001) ++counterLocation[LGL_PREFERRED];
        else if(ratioJacLGL/ratioJacEq >= 1.001) ++counterLocation[EQUI_PREFERRED];
        else                                     ++counterLocation[NO_PREFERRED_LOCATION];
      }
    }
    END_SU2_OMP_FOR

    /*--- Carry out the reduction over the threads. ---*/
    SU2_OMP_CRITICAL
    {
      for(unsigned short i=0; i<=NO_PREFERRED_LOCATION; ++i)
        counterGridLocation[i] += counterLocation[i];
    }
    END_SU2_OMP_CRITICAL

    /*--- When MPI is used, determine the global values of counterGridLocation.
          Only a single thread needs to do this. ---*/
#ifdef HAVE_MPI
    SU2_OMP_SINGLE
    {
      unsigned long tmpCounter[NO_PREFERRED_LOCATION+1];
      for(unsigned short i=0; i<=NO_PREFERRED_LOCATION; ++i)
        tmpCounter[i] = counterGridLocation[i];

      const int count = NO_PREFERRED_LOCATION+1;
      SU2_MPI::Allreduce(tmpCounter, counterGridLocation, count,
                         MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
    }
    END_SU2_OMP_SINGLE
#endif

    /*--- Terminate if elements with no valid location are present, i.e. elements
          for which both node distributions lead to negative Jacobians. ---*/
    if( counterGridLocation[NO_VALID_LOCATION] ) {

      SU2_OMP_SINGLE
      {
        ostringstream message;
        message << "Found " << counterGridLocation[NO_VALID_LOCATION] << " elements with "
                << "negative Jacobians for both the LGL and equidistand node distribution.";
        SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
      }
      END_SU2_OMP_SINGLE
    }

    /*--- Check if there are both elements present for which the LGL distribution must
          be used and elements for which the equidistant distribution must be used.
          Also in that case, terminate. ---*/
    if(counterGridLocation[LGL_ONLY] && counterGridLocation[EQUI_ONLY]) {

      SU2_OMP_SINGLE
      {
        ostringstream message;
        message << "Found " << counterGridLocation[EQUI_ONLY] 
                << " elements with negative Jacobians for LGL distribution"
                << " and " << counterGridLocation[LGL_ONLY]
                << " elements with negative Jacobians for equidistant distribution.";
        SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
      }
      END_SU2_OMP_SINGLE
    }

    /*--- All elements can have a valid mapping to the standard elements.
          Determine the required/preferred node distribution. ---*/
    ENUM_FEM_GRID_LOCATION_CONFIG gridLocation;
    if(     counterGridLocation[LGL_ONLY] > 0)
      gridLocation = LGL;
    else if(counterGridLocation[EQUI_ONLY] > 0)
      gridLocation = EQUIDISTANT;
    else if(counterGridLocation[LGL_PREFERRED] > counterGridLocation[EQUI_PREFERRED])
      gridLocation = LGL;
    else
      gridLocation = EQUIDISTANT;

    /*--- Check for a possible override from the user. ---*/
    if(config->GetKind_FEM_GridDOFsLocation() == LGL) {
      if(counterGridLocation[EQUI_ONLY] > 0) {

        SU2_OMP_SINGLE
        SU2_MPI::Error(string("User specified to use LGL grid DOFs, but only equidistant DOFs are valid"),
                       CURRENT_FUNCTION);
        END_SU2_OMP_SINGLE
      }

      gridLocation = LGL;
    }

    if(config->GetKind_FEM_GridDOFsLocation() == EQUIDISTANT) {
      if(counterGridLocation[LGL_ONLY] > 0) {

        SU2_OMP_SINGLE
        SU2_MPI::Error(string("User specified to use Equidistant grid DOFs, but only LGL DOFs are valid"),
                       CURRENT_FUNCTION);
        END_SU2_OMP_SINGLE
      }

      gridLocation = EQUIDISTANT;
    }

    /*--- Store the grid location to be used in config
          and store this info also in a boolean. ---*/
    config->SetKind_FEM_GridDOFsLocation(gridLocation);
    const bool useLGL = gridLocation == LGL;

    /*--- Allocate the memory for the matrices to compute the face metrics. ---*/
    AllocateMemoryMatricesMetrics(standardFaceElements, nDim-1, &matricesJacobians,
                                  &matricesCoor, &matricesNormals, &matricesDerCoor);

    /*--- Loop over the local volume elements to determine whether or not the Jacobian
          of the element is constant, to determine whether the Jacobians of boundary
          faces are constant and to determine a length scale for the element. ---*/
    SU2_OMP_FOR_DYN(omp_chunk_size)
    for(unsigned long i=0; i<nElem; ++i) {

      /*--- Determine the order of the polynomials that must be integrated exactly.
            As the intention of this function is to find out whether or not the
            Jacobians are constant, it suffices to look at the polynomial degree
            of the grid. ---*/
      const unsigned short nPolyGrid  = elem[i]->GetNPolyGrid();
      const unsigned short orderExact = config->GetOrderExactIntegrationFEM(nPolyGrid, false);

      /*--- Determine the standard element of the volume element. ---*/
      unsigned long ii;
      for(ii=0; ii<standardVolumeElements.size(); ++ii)
        if( standardVolumeElements[ii]->SameStandardElement(elem[i]->GetVTK_Type(),
                                                            nPolyGrid, orderExact) )
          break;

      /*--- Determine the minimum Jacobian and the Jacobian ratio for
            the point distribution used. ---*/
      su2double jacVolMin, ratioJac;
      if(gridLocation == LGL) {
        jacVolMin = jacMinLGL[i];
        ratioJac  = jacMaxLGL[i]/jacMinLGL[i];
      }
      else {
        jacVolMin = jacMinEqui[i];
        ratioJac  = jacMaxEqui[i]/jacMinEqui[i];
      }

      /*--- Determine whether or not the Jacobian can be considered constant. ---*/
      bool constJacobian = (ratioJac <= 1.000001);
      elem[i]->SetJacobianConsideredConstant(constJacobian);

      /*--- Determine the number of faces for this element. ---*/
      const unsigned short nFaces = standardVolumeElements[ii]->GetNFaces();

      /*--- Initialize the array, which stores whether or not the faces are
            considered to have a constant Jacobian. ---*/
      elem[i]->InitializeJacobianConstantFaces(nFaces);

      /*--- Loop over the faces of this element. ---*/
      su2double jacFaceMax = 0.0;
      for(unsigned short j=0; j<nFaces; ++j) {

        /*--- Determine the VTK type of the face element. ---*/
        const unsigned short VTK_Type = standardVolumeElements[ii]->GetVTK_Face(j);

        /*--- Determine the standard element for this face. ---*/
        unsigned long jj;
        for(jj=0; jj<standardFaceElements.size(); ++jj)
          if( standardFaceElements[jj]->SameStandardElement(VTK_Type, nPolyGrid, orderExact) )
            break;

        /*--- Set the pointer to store the local face connectivity of this face. ---*/
        const unsigned short *connFace = standardVolumeElements[ii]->GetGridConnFace(j);

        /*--- Copy the coordinates into matricesCoor[jj], such that the gemm
              routines can be used to compute the derivatives. ---*/
        const unsigned short nDOFs = standardFaceElements[jj]->GetNDOFs();

        for(unsigned short l=0; l<nDOFs; ++l) {
          const unsigned long nodeID = elem[i]->GetNode(connFace[l]);

          map<unsigned long,unsigned long>::const_iterator MI = globalPointIDToLocalInd.find(nodeID);
          const unsigned long ind = MI->second;
          for(unsigned short k=0; k<nDim; ++k)
            matricesCoor[jj](l,k) = nodes->GetCoord(ind, k);
        }

        /*--- Determine the minimum and maximum value of the face jacobian in the integration
              points and the minimum value of the cosine of the angle between the outward normals
              of the integration points. ---*/
        su2double jacMin, jacMax, cosAngleMin;
        standardFaceElements[jj]->MinMaxFaceJacobians(useLGL,  matricesCoor[jj],
                                                      matricesDerCoor[jj],
                                                      matricesNormals[jj],
                                                      matricesJacobians[jj],
                                                      jacMin, jacMax, cosAngleMin);

        /*--- Compute the ratio between the maximum and minimum Jacobian and determine
              whether the face is considered to have a constant Jacobian. ---*/
        const su2double ratioJacMax = jacMax/jacMin;
        constJacobian = cosAngleMin >= 0.999999 && ratioJacMax <= 1.000001;

        elem[i]->SetJacobianConstantFace(constJacobian, j);

        /*--- Update the maximum value of the Jacobian of all faces surrounding
              the current element. ---*/
        jacFaceMax = max(jacFaceMax, jacMax);
      }

      /*--- Determine the length scale of the element. This is needed to compute the
            computational weights of an element when time accurate local time
            stepping is employed and to determine a tolerance when periodic
            transformations are present. Note that a factor 2 must be taken into
            account, which is the length scale of all the reference elements used
            in this code. ---*/
      const su2double lenScale = 2.0*jacVolMin/jacFaceMax;
      elem[i]->SetLengthScale(lenScale);
    }
    END_SU2_OMP_FOR

  }
  END_SU2_OMP_PARALLEL
}

void CPhysicalGeometry::DetermineFEMGraphWeights(CConfig                                    *config,
                                                 const vector<CFaceOfElement>               &localFaces,
                                                 const vector<vector<unsigned long> >       &adjacency,
                                                 const map<unsigned long, CUnsignedShort2T> &mapExternalElemIDToTimeLevel,
                                                 vector<passivedouble>                      &vwgt,
                                                 vector<vector<passivedouble> >             &adjwgt) {

  /*--- Define the linear partitioning of the elements. ---*/
  CLinearPartitioner elemPartitioner(Global_nElem, 0);

  /*--- Determine the maximum time level that occurs in the grid. ---*/
  unsigned short maxTimeLevel = 0;
  for(unsigned long i=0; i<nElem; ++i)
    maxTimeLevel = max(maxTimeLevel, elem[i]->GetTimeLevel());

#ifdef HAVE_MPI
  unsigned short maxTimeLevelLocal = maxTimeLevel;
  SU2_MPI::Allreduce(&maxTimeLevelLocal, &maxTimeLevel, 1,
                     MPI_UNSIGNED_SHORT, MPI_MAX, SU2_MPI::GetComm());
#endif

  /*--- Allocate the memory to store the weights of the graph. ---*/
  vwgt.resize(2*nElem);

  adjwgt.resize(nElem);
  for(unsigned long i=0; i<nElem; ++i)
    adjwgt[i].resize(adjacency[i].size());

  /*--------------------------------------------------------------------------*/
  /* Step 1: Determine the vertex weights of the graph. Per element two       */
  /*         weights are determined. The first weight is proportional to the  */
  /*         amount of work for the volume element. The second weight is the  */
  /*         number of DOFs of the element, such that the number of DOFs per  */
  /*         rank will also be the same.                                      */
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the elements to determine the amount of computational work.
        This amount has a contribution from both the volume integral and
        surface integral to allow for Discontinuous and Continuous Galerkin
        schemes. For the latter the contribution of the surface integral will
        be negligible to the total amount of work. ---*/
  for(unsigned long i=0; i<nElem; ++i) {

    /*--- Easier storage of the two indices for the vertex weights of this element. ---*/
    const unsigned long ind0 = 2*i;
    const unsigned long ind1 = ind0 + 1;

    /*--- Determine whether or not the Jacobian is constant, get the VTK
          type and the polynomial degree of the solution, and determine
          the polynomial degree that must be integrated exactly. ---*/
    bool constJac = elem[i]->GetJacobianConsideredConstant();

    unsigned short VTK_Type   = elem[i]->GetVTK_Type();
    unsigned short nPolySol   = elem[i]->GetNPolySol();
    unsigned short orderExact = config->GetOrderExactIntegrationFEM(nPolySol, constJac);;

    /*--- Determine the corresponding standard volume element. ---*/
    unsigned long ii;
    for(ii=0; ii<standardVolumeElements.size(); ++ii)
      if( standardVolumeElements[ii]->SameStandardElement(VTK_Type, nPolySol, orderExact) )
        break;

    /*--- Initialize the computational work for this element to the work associated
          with the volume. The computational work is stored in the 1st vertex weight. ---*/
    vwgt[ind0] = standardVolumeElements[ii]->WorkEstimateVolume(config);

    /*------------------------------------------------------------------------*/
    /*--- Determine the computational weight of the surface integral in    ---*/
    /*--- the DG-FEM formulation.                                          ---*/
    /*------------------------------------------------------------------------*/

    if(config->GetKind_FEM_Flow() == DG) {

      /*--- Determine the global element ID of this element. ---*/
      const unsigned long elemID = i + elemPartitioner.GetFirstIndexOnRank(rank);

      /*--- Get the global IDs of the corner points of all the faces of this element. ---*/
      unsigned short nFaces;
      unsigned short nPointsPerFace[6];
      unsigned long  faceConn[6][4];

      elem[i]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

      /*--- Loop over the number of faces of this element. ---*/
      for(unsigned short j=0; j<nFaces; ++j) {

        /*--- Determine the VTK type of the face element, whether or not the
              Jacobian of the face is constant and the order of the polynomials
              that must be integrated exactly. ---*/
        VTK_Type   = standardVolumeElements[ii]->GetVTK_Face(j);
        constJac   = elem[i]->GetJacobianConstantFace(j);
        orderExact = config->GetOrderExactIntegrationFEM(nPolySol, constJac);

        /*--- Determine the standard element for this face. ---*/
        unsigned long jj;
        for(jj=0; jj<standardFaceElements.size(); ++jj)
          if( standardFaceElements[jj]->SameStandardElement(VTK_Type, nPolySol, orderExact) )
            break;

        /*--- Create an object of the class CFaceOfElement corresponding to
              the data of this face. ---*/
        CFaceOfElement thisFace;
        thisFace.nCornerPoints = nPointsPerFace[j];
        for(unsigned short k=0; k<nPointsPerFace[j]; ++k)
          thisFace.cornerPoints[k] = faceConn[j][k];
        thisFace.CreateUniqueNumbering();

        /*--- Search for this face in localFaces, which contains the
              internal faces of the grid. ---*/
        vector<CFaceOfElement>::const_iterator low;
        low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);

        bool thisFaceFound = false;
        if(low != localFaces.end()) {
          if( !(thisFace < *low) ) thisFaceFound = true;
        }

        /*--- Check if the face is found. ---*/
        if( thisFaceFound ) {

          /*--- This is an internal face. Determine if it is owned by this element. ---*/
          bool faceIsOwned;
          if(elemID == low->elemID0) faceIsOwned =  low->elem0IsOwner;
          else                       faceIsOwned = !low->elem0IsOwner;

          /*--- If this face is owned, update the workload accordingly. ---*/
          if( faceIsOwned )
            vwgt[ind0] += standardFaceElements[jj]->WorkEstimateInternalFace(config,
                                                                             low->elemType0, low->nPolySol0,
                                                                             low->elemType1, low->nPolySol1);
        }
        else {

          /*--- This is a boundary face, which is owned by definition.
                Update the workload accordingly. ---*/
          vwgt[ind0] += standardFaceElements[jj]->WorkEstimateBoundaryFace(config, elem[i]->GetVTK_Type());
        }
      }
    }

    /*--- Set the value of the second vertex weight to the number of
          solution DOFs of the element. ---*/
    vwgt[ind1] = elem[i]->GetNDOFsSol();
  }

  /*--------------------------------------------------------------------------*/
  /*--- Check for boundary faces for which a wall function treatment       ---*/
  /*--- must be used. This type of boundary condition treatment is         ---*/
  /*--- computationally intensive and must be added to the vertex weight   ---*/
  /*--- of the corresponding element.                                      ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the markers and select the ones for which a wall function
        treatment must be carried out. ---*/
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {

    switch (config->GetMarker_All_KindBC(iMarker)) {
      case ISOTHERMAL:
      case HEAT_FLUX: {
        const string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if(config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::NONE) {

          /*--- Retrieve the integer information for this boundary marker.
                The number of points in normal direction for the wall function
                treatment is the first element of this array. */
          const unsigned short *shortInfo = config->GetWallFunction_IntInfo(Marker_Tag);

          /*--- Loop over the local boundary elements for this marker. ---*/
          for(unsigned long l=0; l<nElem_Bound[iMarker]; ++l) {

            /*--- Retrieve the information from the boundary face, such that
                  the corresponding standard element can be determined. ---*/
            const unsigned short VTK_Type   = bound[iMarker][l]->GetVTK_Type();
            const unsigned long  elemID     = bound[iMarker][l]->GetDomainElement()
                                            - elemPartitioner.GetFirstIndexOnRank(rank);
            const unsigned short nPolySol   = elem[elemID]->GetNPolySol();
            const bool constJac             = bound[iMarker][l]->GetJacobianConsideredConstant();
            const unsigned short orderExact = config->GetOrderExactIntegrationFEM(nPolySol, constJac);

            /*--- Determine the corresponding standard element of the face. ---*/
            unsigned long jj;
            for(jj=0; jj<standardFaceElements.size(); ++jj)
              if( standardFaceElements[jj]->SameStandardElement(VTK_Type, nPolySol, orderExact) )
                break;

            /*-- Update the computational work for the corresponding volume element,
                 i.e. the 1st vertex weight. ---*/
            vwgt[2*elemID] += standardFaceElements[jj]->WorkEstimateWallFunctions(config, shortInfo[0],
                                                                                  elem[elemID]->GetVTK_Type());
          }
        }

        break;
      }

      default:  // Just to avoid a compiler warning.
        break;
    }
  }

  /*--- Take the time level into account for the computational weight.
        Note that this correction is only relevant when time  accurate
        local time stepping is employed. ---*/
  for(unsigned long i=0; i<nElem; ++i) {
    const unsigned short diffLevel = maxTimeLevel - elem[i]->GetTimeLevel();
    vwgt[2*i] *= pow(2, diffLevel);
  }

  /*--------------------------------------------------------------------------*/
  /*--- The final weight is obtained by applying a scaling, such that the  ---*/
  /*--- conversion to integer weights (ParMETIS is using integers for the  ---*/
  /*--- weights) does not lead to a significant increase in the load       ---*/
  /*--- imbalance of the computational work.                               ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Determine the minimum of the workload of the elements, i.e. 1st vertex
        weight, over the entire domain. ---*/
  passivedouble minvwgt = vwgt[0];
  for(unsigned long i=0; i<nElem; ++i) minvwgt = min(minvwgt, vwgt[2*i]);

#ifdef HAVE_MPI
  su2double locminvwgt = minvwgt, globminvwgt;
  SU2_MPI::Allreduce(&locminvwgt, &globminvwgt, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
  minvwgt = SU2_TYPE::GetValue(globminvwgt);
#endif

  /*--- Apply the scaling. ---*/
  minvwgt = 100.0/minvwgt;
  for(unsigned long i=0; i<nElem; ++i) vwgt[2*i] *= minvwgt;

  /*--------------------------------------------------------------------------*/
  /* Step 2: Determine the adjacency weights, which are proportional to the   */
  /*         amount of communication needed when the two neighboring          */
  /*         elements are stored on different ranks.                          */
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the number of elements. ---*/
  for(unsigned long i=0; i<nElem; ++i) {

    /*--- Easier storage of the time level of the current element and the
          number of solution DOFs. ---*/
    const unsigned short timeLevel0 = elem[i]->GetTimeLevel();
    const unsigned short nDOFs0     = elem[i]->GetNDOFsSol();

    /*--- Loop over the number of entries in the graph for this element. ---*/
    for(unsigned long j=0; j<adjacency[i].size(); ++j) {

      /*--- Define the variables for the time level and number of solution
            DOFs for the neighboring element. ---*/
      unsigned short timeLevel1, nDOFs1;

      /*--- Check if the neighboring element is stored locally. ---*/
      if(elemPartitioner.GetRankContainingIndex(adjacency[i][j]) == static_cast<unsigned long>(rank)) {

        /*--- Locally stored element. Determine its local index and retrieve
              the time level and number of DOFs. ---*/
        const unsigned long elemID1 = adjacency[i][j] - elemPartitioner.GetFirstIndexOnRank(rank);
        timeLevel1 = elem[elemID1]->GetTimeLevel();
        nDOFs1     = elem[elemID1]->GetNDOFsSol();
      }
      else {

        /*--- The neighbor is an external element. Find it in mapExternalElemIDToTimeLevel
              and set the time level and number of solution DOFs accordingly. ---*/
        map<unsigned long,CUnsignedShort2T>::const_iterator MI;
        MI = mapExternalElemIDToTimeLevel.find(adjacency[i][j]);
        if(MI == mapExternalElemIDToTimeLevel.end())
          SU2_MPI::Error("Entry not found in mapExternalElemIDToTimeLevel", CURRENT_FUNCTION);
        timeLevel1 = MI->second.short0;
        nDOFs1     = MI->second.short1;
      }

      /*--- Determine the difference of the maximum time level that occurs and
            the minimum of the time level of the current and the adjacent element.
            This value can only be nonzero when time accurate local time stepping
            is employed. ---*/
      const unsigned short diffLevel = maxTimeLevel - min(timeLevel0, timeLevel1);

      /*--- Set the edge weight. As ParMetis expects an undirected graph, set the edge weight
            to the sum of the number of DOFs on both sides, multiplied by the weight factor
            to account for different time levels. ---*/
      adjwgt[i][j] = pow(2, diffLevel)*(nDOFs0 + nDOFs1);
    }
  }
}

void CPhysicalGeometry::DetermineFEMStandardElements(CConfig *config) {

  /*--- Vector of three unsigned shorts per entity to determine the different
        element types in the locally stored volume elements. ---*/
  vector<CUnsignedShort3T> elemTypes(5*nElem);

  /*--- Loop over the local elements to fill the entries of elemTypes. ---*/
  for(unsigned long i=0; i<nElem; ++i) {

    /*--- Retrieve the polynomial degree of the grid and solution. ---*/
    const unsigned short nPolyGrid = elem[i]->GetNPolyGrid();
    const unsigned short nPolySol  = elem[i]->GetNPolySol();

    /*--- Determine the order of the polynomials that must be integrated exactly
          for a) the grid when the Jacobian is not constant, b) the solution
          when the Jacobian is constant and c) the solution when the Jacobian
          is not constant. a) is needed to determine whether or not the Jacobian
          is actually constant and b) and c) are needed to create the face
          elements to determine the integration points for the wall model
          as well as to estimate the amount of work per element. ---*/
    const unsigned short orderExactA = config->GetOrderExactIntegrationFEM(nPolyGrid, false);
    const unsigned short orderExactB = config->GetOrderExactIntegrationFEM(nPolySol,  true);
    const unsigned short orderExactC = config->GetOrderExactIntegrationFEM(nPolySol,  false);

    /*--- Store the five variants in elemTypes. ---*/
    unsigned long ind = 5*i;
    elemTypes[ind].short0 = elem[i]->GetVTK_Type();
    elemTypes[ind].short1 = nPolyGrid;
    elemTypes[ind].short2 = orderExactA;

    ++ind;
    elemTypes[ind].short0 = elem[i]->GetVTK_Type();
    elemTypes[ind].short1 = nPolyGrid;
    elemTypes[ind].short2 = orderExactB;

    ++ind;
    elemTypes[ind].short0 = elem[i]->GetVTK_Type();
    elemTypes[ind].short1 = nPolyGrid;
    elemTypes[ind].short2 = orderExactC;

    ++ind;
    elemTypes[ind].short0 = elem[i]->GetVTK_Type();
    elemTypes[ind].short1 = nPolySol;
    elemTypes[ind].short2 = orderExactB;

    ++ind;
    elemTypes[ind].short0 = elem[i]->GetVTK_Type();
    elemTypes[ind].short1 = nPolySol;
    elemTypes[ind].short2 = orderExactC;
  }

  /*--- Sort elemTypesGrid in increasing order and remove the multiple entities.
        Allocate the memory for standardVolumeElements afterwards. ---*/
  sort(elemTypes.begin(), elemTypes.end());
  vector<CUnsignedShort3T>::iterator lastEntry = unique(elemTypes.begin(), elemTypes.end());
  elemTypes.erase(lastEntry, elemTypes.end());

  standardVolumeElements.resize(elemTypes.size(), nullptr);

  /*--- Loop over the standard volume elements and allocate the appropriate objects. ---*/
  for(unsigned long i=0; i<standardVolumeElements.size(); ++i) {

    /*--- Abbreviate the element type, polynomial degree and polynomial order that
          must be integrated exactly for readability. ---*/
    const unsigned short VTK_Type   = elemTypes[i].short0;
    const unsigned short nPoly      = elemTypes[i].short1;
    const unsigned short orderExact = elemTypes[i].short2;

    /*--- Determine the element type and allocate the appropriate object. ---*/
    switch( VTK_Type ) {
      case TRIANGLE:
        standardVolumeElements[i] = new CFEMStandardTriPartition(nPoly, orderExact, false);
        break;
      case QUADRILATERAL:
        standardVolumeElements[i] = new CFEMStandardQuadPartition(nPoly, orderExact, false);
        break;
      case TETRAHEDRON:
        standardVolumeElements[i] = new CFEMStandardTetPartition(nPoly, orderExact);
        break;
      case PYRAMID:
        standardVolumeElements[i] = new CFEMStandardPyraPartition(nPoly, orderExact);
        break;
      case PRISM:
        standardVolumeElements[i] = new CFEMStandardPrismPartition(nPoly, orderExact);
        break;
      case HEXAHEDRON:
        standardVolumeElements[i] = new CFEMStandardHexPartition(nPoly, orderExact);
        break;
      default:  /*--- To avoid a compiler warning. ---*/
        SU2_MPI::Error(string("Unknown volume element. This should not happen"),
                         CURRENT_FUNCTION);
    }
  }

  /*--- Reset elemTypes. ---*/
  elemTypes.clear();

  /*--- Loop over the standard volume elements to determine the standard faces. ---*/
  for(unsigned long i=0; i<standardVolumeElements.size(); ++i) {

    /*--- Loop over the number of different standard faces for this and
          store the VTK type and polynomial degree in elemTypesGrid. ---*/
    for(unsigned short j=0; j<standardVolumeElements[i]->GetnFaceTypes(); ++j)
      elemTypes.push_back(CUnsignedShort3T(standardVolumeElements[i]->GetVTK_TypeFace(j),
                                           standardVolumeElements[i]->GetPolyDegree(),
                                           standardVolumeElements[i]->GetOrderExact()));
  }

  /*--- Sort elemTypes in increasing order and remove the multiple entries.
        Allocate the memory for standardFaceElements afterwards.  ---*/
  sort(elemTypes.begin(), elemTypes.end());
  lastEntry = unique(elemTypes.begin(), elemTypes.end());
  elemTypes.erase(lastEntry, elemTypes.end());

  standardFaceElements.resize(elemTypes.size(), nullptr);

  /*--- Loop over the standard face elements and allocate the appropriate objects. ---*/
  for(unsigned long i=0; i<standardFaceElements.size(); ++i) {

    /*--- Abbreviate the element type, polynomial degree and polynomial order that
          must be integrated exactly for readability. ---*/
    const unsigned short VTK_Type   = elemTypes[i].short0;
    const unsigned short nPoly      = elemTypes[i].short1;
    const unsigned short orderExact = elemTypes[i].short2;

    /*--- Determine the element type and allocate the appropriate object. ---*/
    switch( VTK_Type ) {
      case LINE:
        standardFaceElements[i] = new CFEMStandardLinePartition(nPoly, orderExact);
        break;
      case TRIANGLE:
        standardFaceElements[i] = new CFEMStandardTriPartition(nPoly, orderExact, true);
        break;
      case QUADRILATERAL:
        standardFaceElements[i] = new CFEMStandardQuadPartition(nPoly, orderExact, true);
        break;
      default:  /*--- To avoid a compiler warning. ---*/
        SU2_MPI::Error(string("Unknown surface element. This should not happen"),
                       CURRENT_FUNCTION);
    }
  }
}

void CPhysicalGeometry::DetermineGraphAdjacency(const vector<CFaceOfElement>   &localMatchingFaces,
                                                vector<vector<unsigned long> > &adjacency) {

  /*--- Define the linear partitioning of the elements. ---*/
  CLinearPartitioner elemPartitioner(Global_nElem, 0);

  /*--- Allocate the memory for the first index of adjacency. ---*/
  adjacency.resize(nElem);

  /*--- Loop over the matching faces to create the adjacency
        coming from internal faces. ---*/
  for(vector<CFaceOfElement>::const_iterator FI =localMatchingFaces.begin();
                                             FI!=localMatchingFaces.end(); ++FI) {

    /*--- Determine the local index of elem0, which is always stored locally,
          and add elemID1 to the adjacency list. ---*/
    const unsigned long elem0 = FI->elemID0 - elemPartitioner.GetFirstIndexOnRank(rank);
    adjacency[elem0].push_back(FI->elemID1);

    /*--- Check if this is not a periodic face and if the second element is
          also a local element. If so, add elemID0 to the adjacency list ---*/
    if(FI->periodicIndex == 0) {
      if(elemPartitioner.GetRankContainingIndex(FI->elemID1) == static_cast<unsigned long>(rank)) {
        const unsigned long elem1 = FI->elemID1 - elemPartitioner.GetFirstIndexOnRank(rank);
        adjacency[elem1].push_back(FI->elemID0);
      }
    }
  }

  /*--- It is possible that some neighbors appear multiple times due to e.g.
        periodic boundary conditions. ParMETIS is not able to deal with this
        situation, hence these multiple entries must be removed. ---*/
  for(unsigned long i=0; i<nElem; ++i) {
    sort(adjacency[i].begin(), adjacency[i].end());
    vector<unsigned long>::iterator lastEntry;
    lastEntry = unique(adjacency[i].begin(), adjacency[i].end());
    adjacency[i].erase(lastEntry, adjacency[i].end());
  }

  /*--- Due to periodic boundary conditions it is also possible that self entries
        are present. ParMETIS is not able to deal with self entries, hence
        they must be removed as well. ---*/
  for(unsigned long i=0; i<nElem; ++i) {
    const unsigned long globalElemID = i + elemPartitioner.GetFirstIndexOnRank(rank);
    unsigned long nEntriesNew = adjacency[i].size();

    for(unsigned long j=0; j<adjacency[i].size(); ++j) {
      if(adjacency[i][j] == globalElemID) {
        adjacency[i][j] = ULONG_MAX;
        --nEntriesNew;
      }
    }

    sort(adjacency[i].begin(), adjacency[i].end());
    adjacency[i].resize(nEntriesNew);
  }

  /*--- Possibly add the connectivities in the graph from the wall function
        treatment. As these connectivities are one-sided, they must be stored
        and communicated later. ---*/
  vector<unsigned long> additionalExternalEntriesGraph;

  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
    for(unsigned long l=0; l<nElem_Bound[iMarker]; ++l) {

      /*--- Get the global and local element ID adjacent to this boundary face. ---*/
      const unsigned long globalElemID = bound[iMarker][l]->GetDomainElement();
      const unsigned long elemID       = globalElemID - elemPartitioner.GetFirstIndexOnRank(rank);

      /*--- Get the number of donor elements for the wall function treatment
            and the pointer to the array which stores this info. ---*/
      const unsigned short nDonors = bound[iMarker][l]->GetNDonorsWallFunctions();
      const unsigned long  *donors = bound[iMarker][l]->GetDonorsWallFunctions();

      /*--- Loop over the number of donors and add the entry in the graph,
            if not already present. ---*/
      for(unsigned short i=0; i<nDonors; ++i) {

        if(donors[i] != globalElemID) {
          if( !binary_search(adjacency[elemID].begin(), adjacency[elemID].end(),
                             donors[i]) ) {

            /*--- Donor not present in the graph for elemID. Add it and sort it
                  afterwards, such that a binary search can be applied later. ---*/
            adjacency[elemID].push_back(donors[i]);
            sort(adjacency[elemID].begin(), adjacency[elemID].end());

            /*--- Check if the donor element is stored locally. ---*/
            if(elemPartitioner.GetRankContainingIndex(donors[i]) == static_cast<unsigned long>(rank)) {

              /*--- Donor is stored locally. Add the entry to the graph
                    and sort it afterwards. ---*/
              const unsigned long localDonorID = donors[i] - elemPartitioner.GetFirstIndexOnRank(rank);
              adjacency[localDonorID].push_back(globalElemID);
              sort(adjacency[localDonorID].begin(), adjacency[localDonorID].end());
            }
            else {

              /*--- Donor is stored externally. Store the graph entry in
                    additionalExternalEntriesGraph. ---*/
              additionalExternalEntriesGraph.push_back(donors[i]);
              additionalExternalEntriesGraph.push_back(globalElemID);
            }
          }
        }
      }
    }
  }

#ifdef HAVE_MPI
  /*--- Create the send buffers with the additional graph data and determine
        to which ranks data is sent. ---*/
  vector<vector<unsigned long> > sendBufsGraphData(size, vector<unsigned long>(0));
  vector<int> sendToRank(size, 0);

  /*--- Loop over the additional entries. Note that both the donor and the
        element adjacent to the wall boundary is stored in this vector. ---*/
  for(unsigned long i=0; i<additionalExternalEntriesGraph.size(); i+=2) {

    /*--- Determine the rank where this external is stored and update
          the corresponding communication buffers accordingly. ---*/
    const unsigned long rankElem = elemPartitioner.GetRankContainingIndex(additionalExternalEntriesGraph[i]);

    sendBufsGraphData[rankElem].push_back(additionalExternalEntriesGraph[i]);
    sendBufsGraphData[rankElem].push_back(additionalExternalEntriesGraph[i+1]);
    sendToRank[rankElem] = 1;
  }

  /*-- Determine to how many ranks this rank will send data and from how
       many ranks it will receive data. ---*/
  int nRankSend = 0;
  for(int i=0; i<size; ++i) {
    if( sendToRank[i] ) ++nRankSend;
  }

  int nRankRecv;
  vector<int> sizeSend(size, 1);
  SU2_MPI::Reduce_scatter(sendToRank.data(), &nRankRecv, sizeSend.data(),
                          MPI_INT, MPI_SUM, SU2_MPI::GetComm());

  /*--- Send the data using non-blocking sends. ---*/
  vector<SU2_MPI::Request> sendReqs(nRankSend);
  nRankSend = 0;
  for(int i=0; i<size; ++i) {
    if( sendToRank[i] )
      SU2_MPI::Isend(sendBufsGraphData[i].data(), sendBufsGraphData[i].size(),
                     MPI_UNSIGNED_LONG, i, i, SU2_MPI::GetComm(),
                     &sendReqs[nRankSend++]);
  }

  /*--- Loop over the number of ranks from which this rank receives data. ---*/
  for(int i=0; i<nRankRecv; ++i) {

    /*--- Block until a message with unsigned longs arrives from any processor.
          Determine the source and the size of the message. ---*/
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    /*--- Allocate the memory for the receive buffer and receive the message
          using a blocking receive. ---*/
    vector<unsigned long> recvBuf(sizeMess);
    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG,
                   source, rank, SU2_MPI::GetComm(), &status);

    /*--- Loop over the contents of the receive buffer and update the
          graph accordingly. ---*/
    for(int j=0; j<sizeMess; j+=2) {
      const unsigned long elemID = recvBuf[j] - elemPartitioner.GetFirstIndexOnRank(rank);
      adjacency[elemID].push_back(recvBuf[j+1]);
      sort(adjacency[elemID].begin(), adjacency[elemID].end());
    }
  }

  /*--- Complete the non-blocking sends amd synchronize the ranks, because
        wild cards have been used in the above communication. ---*/
  SU2_MPI::Waitall(nRankSend, sendReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Barrier(SU2_MPI::GetComm());
#endif
}

void CPhysicalGeometry::DetermineMatchingFacesFEMGrid(const CConfig          *config,
                                                      vector<CFaceOfElement> &localFaces) {

  /*--- Loop over all elements to determine the faces of the elements. ---*/
  for(unsigned long k=0; k<nElem; ++k) {

    /*--- Get the global IDs of the corner points of all the faces of this element. ---*/
    unsigned short nFaces;
    unsigned short nPointsPerFace[6];
    unsigned long  faceConn[6][4];

    elem[k]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

    /*--- Loop over the faces of the element and add them to localFaces. ---*/
    for(unsigned short i=0; i<nFaces; ++i) {

      /*--- Set the information for this face. ---*/
      CFaceOfElement thisFace;
      thisFace.nCornerPoints = nPointsPerFace[i];
      for(unsigned short j=0; j<nPointsPerFace[i]; ++j)
        thisFace.cornerPoints[j] = faceConn[i][j];
      thisFace.elemID0    = elem[k]->GetGlobalElemID();
      thisFace.nPolySol0  = elem[k]->GetNPolySol();
      thisFace.nDOFsElem0 = elem[k]->GetNDOFsSol();
      thisFace.elemType0  = elem[k]->GetVTK_Type();

      /*--- Create a unique number for the face and add it to localFaces. ---*/
      thisFace.CreateUniqueNumbering();
      localFaces.push_back(thisFace);
    }
  }

  /*--- Sort localFaces in increasing order. */
  sort(localFaces.begin(), localFaces.end());

  /*--- Loop over the boundary markers and its faces in order to flag the
        physical boundary faces in localFaces. Note that use is made of the
        overloaded function GetCornerPointsAllFaces, which explains the
        dimensions of the variables used in the function call. Also note that
        the periodic boundaries are excluded, because they are not physical.
        The face is invalidated by setting the element ID to an invalid value. ---*/
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
    if(config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY) {
      for(unsigned long k=0; k<nElem_Bound[iMarker]; ++k) {

        unsigned short nFaces;
        unsigned short nPointsPerFace[6];
        unsigned long  faceConn[6][4];
        bound[iMarker][k]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

        CFaceOfElement thisFace;
        thisFace.nCornerPoints = nPointsPerFace[0];
        for(unsigned short j=0; j<nPointsPerFace[0]; ++j)
          thisFace.cornerPoints[j] = faceConn[0][j];
        thisFace.CreateUniqueNumbering();

        vector<CFaceOfElement>::iterator low;
        low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);

        low->elemID0 = Global_nElem + 10;
      }
    }
  }

  /*--- Loop over the faces and check for double entries. If a double entry is
        found, the elemID from the second entry is copied to the first entry,
        the polynomial degree is updated, and the second entry is invalidated. ---*/
  unsigned long nFacesLoc = localFaces.size();
  for(unsigned long i=1; i<nFacesLoc; ++i) {
    if(localFaces[i] == localFaces[i-1]) {
      localFaces[i-1].elemID1    = localFaces[i].elemID0;
      localFaces[i-1].nPolySol1  = localFaces[i].nPolySol0;
      localFaces[i-1].nDOFsElem1 = localFaces[i].nDOFsElem0;
      localFaces[i-1].elemType1  = localFaces[i].elemType0;
      localFaces[i].elemID0      = Global_nElem + 10;
    }
  }

  /*--- Remove the invalidated faces. This is accomplished by giving the
        face four points a global node ID that is larger than the largest
        point ID in the grid. In this way the sorting operator puts
        these faces at the end of the vector, see also the < operator
        of CFaceOfElement. ---*/
  unsigned long nFacesLocOr = nFacesLoc;
  for(unsigned long i=0; i<nFacesLocOr; ++i) {
    if(localFaces[i].elemID0 > Global_nElem) {
      localFaces[i].nCornerPoints = 4;
      localFaces[i].cornerPoints[0] = Global_nPoint;
      localFaces[i].cornerPoints[1] = Global_nPoint;
      localFaces[i].cornerPoints[2] = Global_nPoint;
      localFaces[i].cornerPoints[3] = Global_nPoint;
      --nFacesLoc;
    }
  }

  sort(localFaces.begin(), localFaces.end());
  localFaces.resize(nFacesLoc);

  /*--- The remaining part of this function is only needed if MPI is used. ---*/ 
#ifdef HAVE_MPI

  /*--- Determine and store the faces, for which only one neighboring element
        was found. For these faces the other neighbor might be stored on a
        different rank (unless there are periodic or non-matching interfaces). ---*/  
  vector<CFaceOfElement> localFacesComm;
  for(unsigned long i=0; i<nFacesLoc; ++i)
    if(localFaces[i].elemID1 > Global_nElem) localFacesComm.push_back(localFaces[i]);

  /*--- Determine the maximum global point ID that occurs in localFacesComm
        of all ranks. Note that only the first point is taken into account,
        because this point determines the rank where the face is sent to. ---*/
  unsigned long nFacesLocComm = localFacesComm.size();
  unsigned long maxPointIDLoc = 0;
  for(unsigned long i=0; i<nFacesLocComm; ++i)
    maxPointIDLoc = max(maxPointIDLoc, localFacesComm[i].cornerPoints[0]);

  unsigned long maxPointID;
  SU2_MPI::Allreduce(&maxPointIDLoc, &maxPointID, 1, MPI_UNSIGNED_LONG,
                     MPI_MAX, SU2_MPI::GetComm());
  ++maxPointID;

  /*--- Create a linear partition for the points that occur in
        the faces of localFacesComm. ---*/
  CLinearPartitioner facePartitioner(maxPointID,0);  

  /*--- Define the send buffers for the faces. ---*/
  vector<vector<unsigned long> > sendBufFace(size, vector<unsigned long>(0));

  /*--- Loop over the local faces to be communicated and fill
        the entries of sendBufFace. ---*/
  for(unsigned long i=0; i<nFacesLocComm; ++i) {
    const unsigned long rankFace = facePartitioner.GetRankContainingIndex(localFacesComm[i].cornerPoints[0]);

    sendBufFace[rankFace].push_back(localFacesComm[i].nCornerPoints);
    sendBufFace[rankFace].push_back(localFacesComm[i].cornerPoints[0]);
    sendBufFace[rankFace].push_back(localFacesComm[i].cornerPoints[1]);
    sendBufFace[rankFace].push_back(localFacesComm[i].cornerPoints[2]);
    sendBufFace[rankFace].push_back(localFacesComm[i].cornerPoints[3]);
    sendBufFace[rankFace].push_back(localFacesComm[i].elemID0);
    sendBufFace[rankFace].push_back(localFacesComm[i].nPolySol0);
    sendBufFace[rankFace].push_back(localFacesComm[i].nDOFsElem0);
    sendBufFace[rankFace].push_back(localFacesComm[i].elemType0);
  }

  /*--- Determine the number of ranks from which I receive a message. */
  vector<unsigned long> counter(size, 0);
  unsigned long nMessSend = 0;
  for(int i=0; i<size; ++i) {
    if( sendBufFace[i].size() ) {counter[i] = 1; ++nMessSend;}
  }

  vector<int> sizeRecv(size, 1);

  unsigned long nMessRecv;
  SU2_MPI::Reduce_scatter(counter.data(), &nMessRecv, sizeRecv.data(),
                          MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  /*--- Send the data using nonblocking sends. ---*/
  vector<SU2_MPI::Request> commReqs(max(nMessSend,nMessRecv));

  nMessSend = 0;
  for(int i=0; i<size; ++i) {
    if( sendBufFace[i].size() ) {
      SU2_MPI::Isend(sendBufFace[i].data(), sendBufFace[i].size(),
                     MPI_UNSIGNED_LONG, i, i, SU2_MPI::GetComm(),
                     &commReqs[nMessSend]);
      ++nMessSend;
    }
  }

  /*--- Loop over the number of ranks from which faces are received.
        Receive the messages and store the faces in facesRecv. ---*/
  vector<vector<CFaceOfElement> > facesRecv(nMessRecv, vector<CFaceOfElement>(0));
  vector<int> rankRecv(nMessRecv);

  for(unsigned long i=0; i<nMessRecv; ++i) {
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    rankRecv[i] = status.MPI_SOURCE;
    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    vector<unsigned long> recvBuf(sizeMess);
    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG,
                  rankRecv[i], rank, SU2_MPI::GetComm(), &status);

    const int nFacesRecv = sizeMess/9;
    facesRecv[i].resize(nFacesRecv);

    unsigned long ii = 0;
    for(int j=0; j<nFacesRecv; ++j, ii+=9) {
      facesRecv[i][j].nCornerPoints   = recvBuf[ii];
      facesRecv[i][j].cornerPoints[0] = recvBuf[ii+1];
      facesRecv[i][j].cornerPoints[1] = recvBuf[ii+2];
      facesRecv[i][j].cornerPoints[2] = recvBuf[ii+3];
      facesRecv[i][j].cornerPoints[3] = recvBuf[ii+4];
      facesRecv[i][j].elemID0         = recvBuf[ii+5];
      facesRecv[i][j].nPolySol0       = recvBuf[ii+6];
      facesRecv[i][j].nDOFsElem0      = recvBuf[ii+7];
      facesRecv[i][j].elemType0       = recvBuf[ii+8];
    }
  }

  /*--- The exact contents of facesRecv is needed when communicating back
        the data to the original ranks, so the data is copied to
        localFacesComm, which is not needed anymore. ---*/
  localFacesComm.clear();
  for(unsigned long i=0; i<nMessRecv; ++i)
    localFacesComm.insert(localFacesComm.end(), facesRecv[i].begin(), facesRecv[i].end());

  /*--- Sort localFacesComm in increasing order and determine the
        double entries. ---*/
  sort(localFacesComm.begin(), localFacesComm.end());

  nFacesLocComm = localFacesComm.size();
  for(unsigned long i=1; i<localFacesComm.size(); ++i) {
    if(localFacesComm[i] == localFacesComm[i-1]) {
      localFacesComm[i-1].elemID1    = localFacesComm[i].elemID0;
      localFacesComm[i-1].nPolySol1  = localFacesComm[i].nPolySol0;
      localFacesComm[i-1].nDOFsElem1 = localFacesComm[i].nDOFsElem0;
      localFacesComm[i-1].elemType1  = localFacesComm[i].elemType0;

      localFacesComm[i].nCornerPoints = 4;
      localFacesComm[i].cornerPoints[0] = Global_nPoint;
      localFacesComm[i].cornerPoints[1] = Global_nPoint;
      localFacesComm[i].cornerPoints[2] = Global_nPoint;
      localFacesComm[i].cornerPoints[3] = Global_nPoint;
      --nFacesLocComm;
    }
  }

  /*--- Remove the invalidated faces. ---*/
  sort(localFacesComm.begin(), localFacesComm.end());
  localFacesComm.resize(nFacesLocComm);

  /*--- Complete the first round of non-blocking sends. ---*/
  SU2_MPI::Waitall(nMessSend, commReqs.data(), MPI_STATUSES_IGNORE);

  /*--- Send the data back to the requesting ranks. ---*/
  for(unsigned long i=0; i<nMessRecv; ++i) {
    sendBufFace[i].resize(9*facesRecv[i].size());

    unsigned long ii = 0;
    for(unsigned long j=0; j<facesRecv[i].size(); ++j, ii+=9) {

      /*--- Copy the first 5 elements of the received face in the
            send buffer. ---*/
      sendBufFace[i][ii]   = facesRecv[i][j].nCornerPoints;
      sendBufFace[i][ii+1] = facesRecv[i][j].cornerPoints[0];
      sendBufFace[i][ii+2] = facesRecv[i][j].cornerPoints[1];
      sendBufFace[i][ii+3] = facesRecv[i][j].cornerPoints[2];
      sendBufFace[i][ii+4] = facesRecv[i][j].cornerPoints[3];

      /*--- Search for this face in localFacesComm. ---*/
      vector<CFaceOfElement>::iterator low;
      low = lower_bound(localFacesComm.begin(), localFacesComm.end(),
                        facesRecv[i][j]);

      /*--- Fill the remainder of the send buffer, depending on
            whether the received face corresponds to the first
            or second element. ---*/
      if(facesRecv[i][j].elemID0 == low->elemID0) {
        sendBufFace[i][ii+5] = low->elemID1;
        sendBufFace[i][ii+6] = low->nPolySol1;
        sendBufFace[i][ii+7] = low->nDOFsElem1;
        sendBufFace[i][ii+8] = low->elemType1;
      }
      else {
        sendBufFace[i][ii+5] = low->elemID0;
        sendBufFace[i][ii+6] = low->nPolySol0;
        sendBufFace[i][ii+7] = low->nDOFsElem0;
        sendBufFace[i][ii+8] = low->elemType0;
      }
    }

    /*--- Send the data using a non-blocking communication. ---*/
    SU2_MPI::Isend(sendBufFace[i].data(), sendBufFace[i].size(), MPI_UNSIGNED_LONG,
                   rankRecv[i], rankRecv[i]+1, SU2_MPI::GetComm(), &commReqs[i]);
  }

  /*--- Loop over the ranks to which I originally sent my face data.
        The return data contains information about the neighboring element. ---*/
  for(unsigned long i=0; i<nMessSend; ++i) {

    /*--- Wait until a message arrives and determine the source and
          size of the message. ---*/
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank+1, SU2_MPI::GetComm(), &status);
    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);

    /*--- Allocate the memory for the receive buffer and receive the data. ---*/
    vector<unsigned long> recvBuf(sizeMess);
    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG,
                  status.MPI_SOURCE, rank+1, SU2_MPI::GetComm(), &status);

    /*--- Loop over the number of received faces. ---*/
    sizeMess /= 9;
    unsigned long jj = 0;
    for(int j=0; j<sizeMess; ++j, jj+=9) {

      /*--- Create an object of type CFaceOfElement to search in the
            locally stored faces. ---*/
      CFaceOfElement thisFace;
      thisFace.nCornerPoints   = recvBuf[jj];
      thisFace.cornerPoints[0] = recvBuf[jj+1];
      thisFace.cornerPoints[1] = recvBuf[jj+2];
      thisFace.cornerPoints[2] = recvBuf[jj+3];
      thisFace.cornerPoints[3] = recvBuf[jj+4];

      /*-- Search in localFaces for this face and set the information
           of the neighboring element. ---*/
      vector<CFaceOfElement>::iterator low;
      low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);
      low->elemID1    = recvBuf[jj+5];
      low->nPolySol1  = recvBuf[jj+6];
      low->nDOFsElem1 = recvBuf[jj+7];
      low->elemType1  = recvBuf[jj+8];
    }
  }

  SU2_MPI::Waitall(nMessRecv, commReqs.data(), MPI_STATUSES_IGNORE);

  /*--- Wild cards have been used in the communication, so
        synchronize the ranks to avoid problems.          ---*/
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#endif
}

void CPhysicalGeometry::DetermineNonMatchingFacesFEMGrid(const CConfig          *config,
                                                         vector<CFaceOfElement> &localMatchingFaces) {

  /*--- If there are single faces still present in localMatchingFaces, these are
        non-matching faces. Store these faces in singleFaces and remove them
        from localMatchingFaces. ---*/
  vector<CFaceOfElement> singleFaces;

  unsigned long nFacesLoc = localMatchingFaces.size();
  for(unsigned long i=0; i<localMatchingFaces.size(); ++i) {
    if(localMatchingFaces[i].elemID1 > Global_nElem) {

      singleFaces.push_back(localMatchingFaces[i]);

      localMatchingFaces[i].nCornerPoints = 4;
      localMatchingFaces[i].cornerPoints[0] = Global_nPoint;
      localMatchingFaces[i].cornerPoints[1] = Global_nPoint;
      localMatchingFaces[i].cornerPoints[2] = Global_nPoint;
      localMatchingFaces[i].cornerPoints[3] = Global_nPoint;
      --nFacesLoc;
    }
  }

  if( singleFaces.size() ) {
    sort(localMatchingFaces.begin(), localMatchingFaces.end());
    localMatchingFaces.resize(nFacesLoc);
  }

  /*--- Determine the global number of non-matching faces. ---*/
  nFacesLoc = singleFaces.size();
  unsigned long nNonMatchingFaces = nFacesLoc;


#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nFacesLoc, &nNonMatchingFaces, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
#endif

  /*--- If there are no non-matching faces, return, such that in the rest of this
        function it can be assumed that non-matching faces are present. ---*/
  if(nNonMatchingFaces == 0) return;


  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

void CPhysicalGeometry::DetermineOwnershipInternalFaces(vector<CFaceOfElement>               &localFaces,
                                                        map<unsigned long, CUnsignedShort2T> &mapExternalElemIDToTimeLevel) {

  /*--- Define the linear partitioning of the elements. ---*/
  CLinearPartitioner elemPartitioner(Global_nElem, 0);

  /*--- Loop over the locally stored faces. ---*/
  for(vector<CFaceOfElement>::iterator FI =localFaces.begin();
                                       FI!=localFaces.end(); ++FI) {

    /*--- Check for a matching face. ---*/
    if(FI->elemID1 < Global_nElem) {

      /*--- Determine the time level of Elem0, which is always owned. ---*/
      const unsigned long  elemID0    = FI->elemID0 - elemPartitioner.GetFirstIndexOnRank(rank);
      const unsigned short timeLevel0 = elem[elemID0]->GetTimeLevel();

      /*--- Determine the time level of Elem1, which is either owned or
            external. Hence a distinction must be made. ---*/
      unsigned short timeLevel1;
      if(elemPartitioner.GetRankContainingIndex(FI->elemID1) == static_cast<unsigned long>(rank)) {

        const unsigned long elemID1 = FI->elemID1 - elemPartitioner.GetFirstIndexOnRank(rank);
        timeLevel1 = elem[elemID1]->GetTimeLevel();
      }
      else {

        map<unsigned long,CUnsignedShort2T>::const_iterator MI;
        MI = mapExternalElemIDToTimeLevel.find(FI->elemID1);
        if(MI == mapExternalElemIDToTimeLevel.end())
          SU2_MPI::Error("Entry not found in mapExternalElemIDToTimeLevel", CURRENT_FUNCTION);
        timeLevel1 = MI->second.short0;
      }

      /*--- Check if both elements have the same time level. ---*/
      if(timeLevel0 == timeLevel1) {

        /*--- Same time level, hence both elements can own the face. First check whether
              elemID0 == elemID1 (which happens for periodic problems with only one
              element in the periodic direction), because this is a special case. */
        if(FI->elemID0 == FI->elemID1) {

          /*--- This face occurs twice, but should be owned only once. Base this
                decision on the periodic index. ---*/
          FI->elem0IsOwner = FI->periodicIndex < FI->periodicIndexDonor;
        }
        else {

          /*--- Different elements on both sides of the face. The ownership decision
                below makes an attempt to spread the workload evenly. ---*/
          const unsigned long sumElemID = FI->elemID0 + FI->elemID1;
          if( sumElemID%2 )
            FI->elem0IsOwner = FI->elemID0 < FI->elemID1;
          else
            FI->elem0IsOwner = FI->elemID0 > FI->elemID1;
        }
      }
      else {

        /*--- The time level of both elements differ. The element with the smallest
              time level must be the owner of the element. ---*/
        FI->elem0IsOwner = timeLevel0 < timeLevel1;
      }
    }
    else {

      /*--- Non-matching face. Give the ownership to element 0. ---*/
      FI->elem0IsOwner = true;
    }
  }
}

void CPhysicalGeometry::DeterminePeriodicFacesFEMGrid(const CConfig          *config,
                                                      vector<CFaceOfElement> &localFaces) {

  /*--- Determine a mapping from the global point ID to the local index
        of the points.            ---*/
  map<unsigned long,unsigned long> globalPointIDToLocalInd;
  for(unsigned i=0; i<nPoint; ++i) {
    globalPointIDToLocalInd[nodes->GetGlobalIndex(i)] = i;
  }

  /*--- Loop over the number of markers present in the grid and check for a periodic one. ---*/
  for(unsigned short iMarker=0; iMarker<config->GetnMarker_All(); ++iMarker) {
    if(config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {

      /*--- Determine the donor marker and the transformation from the
            current marker to the donor marker.    ---*/
      const unsigned short jMarker = config->GetMarker_Periodic_Donor(config->GetMarker_All_TagBound(iMarker));

      auto center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
      auto angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
      auto trans  = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));

      /*--- Store (center+trans) as it is constant and will be added on. ---*/
      const su2double translation[] = {center[0] + trans[0],
                                       center[1] + trans[1],
                                       center[2] + trans[2]};

      /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
      const su2double theta = angles[0];
      const su2double phi   = angles[1];
      const su2double psi   = angles[2];

      const su2double cosTheta = cos(theta), cosPhi = cos(phi), cosPsi = cos(psi);
      const su2double sinTheta = sin(theta), sinPhi = sin(phi), sinPsi = sin(psi);

      /*--- Compute the rotation matrix. Note that the implicit
            ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
      su2double rotMatrix[3][3];
      rotMatrix[0][0] =  cosPhi*cosPsi;
      rotMatrix[1][0] =  cosPhi*sinPsi;
      rotMatrix[2][0] = -sinPhi;

      rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
      rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
      rotMatrix[2][1] = sinTheta*cosPhi;

      rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
      rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
      rotMatrix[2][2] = cosTheta*cosPhi;

      /*--- Define the vector to store the faces of the donor. Initialize its
            size to the number of local donor faces.    ---*/
      vector<CMatchingFace> facesDonor(nElem_Bound[jMarker]);

      /*------------------------------------------------------------------*/
      /*--- Step 1: Store the information of the local faces of the donor
                    marker in the variables defined above.             ---*/
      /*------------------------------------------------------------------*/

      /*--- Loop over the local elements of the donor marker. ---*/
      for(unsigned long k=0; k<nElem_Bound[jMarker]; ++k) {

        /*--- Get the connectivity of this face. The reason for the used
              function arguments for GetCornerPointsAllFaces, is that
              GetCornerPointsAllFaces is an overloaded function.   ---*/
        unsigned short nFaces;
        unsigned short nPointsPerFace[6];
        unsigned long  faceConn[6][4];
        bound[jMarker][k]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

        /*--- Search for this face in localFaces. It must be present. ---*/
        CFaceOfElement thisFace;
        thisFace.nCornerPoints = nPointsPerFace[0];
        for(unsigned short j=0; j<nPointsPerFace[0]; ++j)
          thisFace.cornerPoints[j] = faceConn[0][j];
        thisFace.CreateUniqueNumbering();

        vector<CFaceOfElement>::iterator low;
        low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);

        /*--- Store the relevant data in facesDonor. ---*/
        facesDonor[k].nDim          = nDim;
        facesDonor[k].nCornerPoints = nPointsPerFace[0];
        facesDonor[k].elemID        = low->elemID0;
        facesDonor[k].nPoly         = low->nPolySol0;
        facesDonor[k].nDOFsElem     = low->nDOFsElem0;
        facesDonor[k].elemType      = low->elemType0;

        for(unsigned short j=0; j<nPointsPerFace[0]; ++j) {
          map<unsigned long,unsigned long>::const_iterator MI;
          MI = globalPointIDToLocalInd.find(faceConn[0][j]);
          unsigned long ind = MI->second;

          for(unsigned l=0; l<nDim; ++l)
            facesDonor[k].cornerCoor[j][l] = nodes->GetCoord(ind, l);
        }

        /*--- Create the tolerance for this face and sort the coordinates. ---*/
        facesDonor[k].SortFaceCoordinates();
      }

      /*------------------------------------------------------------------*/
      /*--- Step 2: In parallel mode the data of the donor marker is
                    gathered on all ranks.                             ---*/
      /*------------------------------------------------------------------*/

#ifdef HAVE_MPI

      /*--- Check if this is indeed a parallel simulation. ---*/
      if(size > 1) {

        /*--- Allocate the memory for the size arrays in Allgatherv. ---*/
        vector<int> recvCounts(size), displs(size);

        /*--- Create the values of recvCounts for the gather of the facesDonor. ---*/
        int sizeLocal = facesDonor.size();

        SU2_MPI::Allgather(&sizeLocal, 1, MPI_INT, recvCounts.data(), 1,
                           MPI_INT, SU2_MPI::GetComm());

        /*--- Create the data for the vector displs from the known values of
              recvCounts. Also determine the total size of the data.   ---*/
        displs[0] = 0;
        for(int i=1; i<size; ++i) displs[i] = displs[i-1] + recvCounts[i-1];

        int sizeGlobal = displs.back() + recvCounts.back();

        /*--- SU2_MPI does not support the communication of derived data types,
              at least not at the moment when this was coded. Therefore the data
              from facesDonor is copied into buffers of primitive types, which
              are communicated. ---*/
        vector<unsigned short> shortLocBuf(5*sizeLocal);
        vector<unsigned long>  longLocBuf(sizeLocal);
        vector<su2double>      doubleLocBuf(13*sizeLocal);

        unsigned long cS=0, cL=0, cD=0;
        for(vector<CMatchingFace>::const_iterator fIt =facesDonor.begin();
                                                  fIt!=facesDonor.end(); ++fIt) {
          shortLocBuf[cS++] = fIt->nCornerPoints;
          shortLocBuf[cS++] = fIt->nDim;
          shortLocBuf[cS++] = fIt->nPoly;
          shortLocBuf[cS++] = fIt->nDOFsElem;
          shortLocBuf[cS++] = fIt->elemType;

          longLocBuf[cL++] = fIt->elemID;

          doubleLocBuf[cD++] = fIt->cornerCoor[0][0];
          doubleLocBuf[cD++] = fIt->cornerCoor[0][1];
          doubleLocBuf[cD++] = fIt->cornerCoor[0][2];

          doubleLocBuf[cD++] = fIt->cornerCoor[1][0];
          doubleLocBuf[cD++] = fIt->cornerCoor[1][1];
          doubleLocBuf[cD++] = fIt->cornerCoor[1][2];

          doubleLocBuf[cD++] = fIt->cornerCoor[2][0];
          doubleLocBuf[cD++] = fIt->cornerCoor[2][1];
          doubleLocBuf[cD++] = fIt->cornerCoor[2][2];

          doubleLocBuf[cD++] = fIt->cornerCoor[3][0];
          doubleLocBuf[cD++] = fIt->cornerCoor[3][1];
          doubleLocBuf[cD++] = fIt->cornerCoor[3][2];

          doubleLocBuf[cD++] = fIt->tolForMatching;
        }

        /*--- Gather the faces from all ranks to all ranks. Use Allgatherv
              for this purpose.                  ---*/
        vector<unsigned short> shortGlobBuf(5*sizeGlobal);
        vector<unsigned long>  longGlobBuf(sizeGlobal);
        vector<su2double>      doubleGlobBuf(13*sizeGlobal);

        SU2_MPI::Allgatherv(longLocBuf.data(), longLocBuf.size(), MPI_UNSIGNED_LONG,
                            longGlobBuf.data(), recvCounts.data(), displs.data(),
                            MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

        for(int i=0; i<size; ++i) {
          recvCounts[i] *= 5; displs[i] *= 5;
        }

        SU2_MPI::Allgatherv(shortLocBuf.data(), shortLocBuf.size(), MPI_UNSIGNED_SHORT,
                            shortGlobBuf.data(), recvCounts.data(), displs.data(),
                            MPI_UNSIGNED_SHORT, SU2_MPI::GetComm());

        for(int i=0; i<size; ++i) {
          recvCounts[i] /=  5; displs[i] /=  5;
          recvCounts[i] *= 13; displs[i] *= 13;
        }

        SU2_MPI::Allgatherv(doubleLocBuf.data(), doubleLocBuf.size(), MPI_DOUBLE,
                            doubleGlobBuf.data(), recvCounts.data(), displs.data(),
                            MPI_DOUBLE, SU2_MPI::GetComm());

        /*--- Copy the data back into facesDonor, which will contain the
              global information after the copies. ---*/
        facesDonor.resize(sizeGlobal);
        cS = cL = cD = 0;

        for(vector<CMatchingFace>::iterator fIt =facesDonor.begin();
                                            fIt!=facesDonor.end(); ++fIt) {
          fIt->nCornerPoints = shortGlobBuf[cS++];
          fIt->nDim          = shortGlobBuf[cS++];
          fIt->nPoly         = shortGlobBuf[cS++];
          fIt->nDOFsElem     = shortGlobBuf[cS++];
          fIt->elemType      = shortGlobBuf[cS++];

          fIt->elemID = longGlobBuf[cL++];

          fIt->cornerCoor[0][0] = doubleGlobBuf[cD++];
          fIt->cornerCoor[0][1] = doubleGlobBuf[cD++];
          fIt->cornerCoor[0][2] = doubleGlobBuf[cD++];

          fIt->cornerCoor[1][0] = doubleGlobBuf[cD++];
          fIt->cornerCoor[1][1] = doubleGlobBuf[cD++];
          fIt->cornerCoor[1][2] = doubleGlobBuf[cD++];

          fIt->cornerCoor[2][0] = doubleGlobBuf[cD++];
          fIt->cornerCoor[2][1] = doubleGlobBuf[cD++];
          fIt->cornerCoor[2][2] = doubleGlobBuf[cD++];

          fIt->cornerCoor[3][0] = doubleGlobBuf[cD++];
          fIt->cornerCoor[3][1] = doubleGlobBuf[cD++];
          fIt->cornerCoor[3][2] = doubleGlobBuf[cD++];

          fIt->tolForMatching = doubleGlobBuf[cD++];
        }
      }
#endif

      /*--- Sort facesDonor in increasing order. ---*/
      sort(facesDonor.begin(), facesDonor.end());

      /*------------------------------------------------------------------*/
      /*--- Step 3: Find for the current marker the required data in the
                    variables for the donor marker.                    ---*/
      /*------------------------------------------------------------------*/

      /*--- Loop over the local faces of this boundary marker. ---*/
      for(unsigned long k=0; k<nElem_Bound[iMarker]; ++k) {

        /*--- Get the connectivity of this face. The reason for the used
              function arguments for GetCornerPointsAllFaces, is that
              GetCornerPointsAllFaces is an overloaded function.   ---*/
        unsigned short nFaces;
        unsigned short nPointsPerFace[6];
        unsigned long  faceConn[6][4];
        bound[iMarker][k]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

        /*--- Search for this face in localFaces. It must be present. ---*/
        CFaceOfElement thisFace;
        thisFace.nCornerPoints = nPointsPerFace[0];
        for(unsigned short j=0; j<nPointsPerFace[0]; ++j)
          thisFace.cornerPoints[j] = faceConn[0][j];
        thisFace.CreateUniqueNumbering();

        vector<CFaceOfElement>::iterator low;
        low = lower_bound(localFaces.begin(), localFaces.end(), thisFace);

        /*--- Indicate that this face is a periodic face. This is accomplished
              by setting periodicIndex to iMarker + 1.       ---*/
        low->periodicIndex = iMarker + 1;

        /*--- Store the data for this boundary element also in thisMatchingFace,
              such that a search can be carried out in donorFaces. Note that the
              periodic transformation must be applied to the coordinates. ---*/
        CMatchingFace thisMatchingFace;
        thisMatchingFace.nDim          = nDim;
        thisMatchingFace.nCornerPoints = nPointsPerFace[0];
        thisMatchingFace.elemID        = low->elemID0;
        thisMatchingFace.nPoly         = low->nPolySol0;
        thisMatchingFace.nDOFsElem     = low->nDOFsElem0;
        thisMatchingFace.elemType      = low->elemType0;

        for(unsigned short j=0; j<nPointsPerFace[0]; ++j) {
          map<unsigned long,unsigned long>::const_iterator MI;
          MI = globalPointIDToLocalInd.find(faceConn[0][j]);
          unsigned long ind = MI->second;
          const su2double *coor = nodes->GetCoord(ind);

          const su2double dx =             coor[0] - center[0];
          const su2double dy =             coor[1] - center[1];
          const su2double dz = nDim == 3 ? coor[2] - center[2] : su2double(0.0);

          thisMatchingFace.cornerCoor[j][0] = rotMatrix[0][0]*dx + rotMatrix[0][1]*dy
                                            + rotMatrix[0][2]*dz + translation[0];
          thisMatchingFace.cornerCoor[j][1] = rotMatrix[1][0]*dx + rotMatrix[1][1]*dy
                                            + rotMatrix[1][2]*dz + translation[1];
          thisMatchingFace.cornerCoor[j][2] = rotMatrix[2][0]*dx + rotMatrix[2][1]*dy
                                            + rotMatrix[2][2]*dz + translation[2];
        }

        /*--- Create the tolerance for this face and sort the coordinates. ---*/
        thisMatchingFace.SortFaceCoordinates();

        /*--- Check if thisMatchingFace is present in facesDonor. If so, set the
              missing information in the face corresponding to the iterator low. ---*/
        vector<CMatchingFace>::const_iterator donorLow;
        donorLow = lower_bound(facesDonor.begin(), facesDonor.end(), thisMatchingFace);

        if(donorLow != facesDonor.end()) {
          if( !(thisMatchingFace < *donorLow) ) {
            low->elemID1            = donorLow->elemID;
            low->nPolySol1          = donorLow->nPoly;
            low->nDOFsElem1         = donorLow->nDOFsElem;
            low->elemType1          = donorLow->elemType;
            low->periodicIndexDonor = jMarker + 1;
          }
        }
      }
    }
  } 
}

void CPhysicalGeometry::DetermineTimeLevelElements(CConfig                              *config,
                                                   const vector<CFaceOfElement>         &localFaces,
                                                   map<unsigned long, CUnsignedShort2T> &mapExternalElemIDToTimeLevel) {

  /*--- Define the linear partitioning of the elements. ---*/
  CLinearPartitioner elemPartitioner(Global_nElem, 0);

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Initialize the map mapExternalElemIDToTimeLevel.           ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Initialize the time level of external elements to zero.
        First the externals from the faces. ---*/
  for(vector<CFaceOfElement>::const_iterator FI =localFaces.begin();
                                             FI!=localFaces.end(); ++FI) {
    if(FI->elemID1 < Global_nElem) {     /*--- Safeguard against non-matching faces. ---*/

      /*--- Check for external element. This is done by checking the
            local elements and if it is not found, it is an external. ---*/
      unordered_map<unsigned long,unsigned long>::const_iterator UMI;
      UMI = Global_to_Local_Elem.find(FI->elemID1);

      if(UMI == Global_to_Local_Elem.end()) {

        /*--- This element is an external element. Store it in the map
              mapExternalElemIDToTimeLevel if not already done so. ---*/
        map<unsigned long,CUnsignedShort2T>::iterator MI;
        MI = mapExternalElemIDToTimeLevel.find(FI->elemID1);
        if(MI == mapExternalElemIDToTimeLevel.end())
          mapExternalElemIDToTimeLevel[FI->elemID1] = CUnsignedShort2T(0,0);
      }
    }
  }

  /*--- Define the communication buffers to send additional externals
        to other ranks. Also define the buffer recvFromRank which is used
        in Reduce_scatter later on. It indicates whether or not a message
        is sent to a certain rank. ---*/
  vector<vector<unsigned long> > sendBufAddExternals(size, vector<unsigned long>(0));
  vector<int> recvFromRank(size, 0);

  /*--- Add the externals from the wall function donors. Loop over the
        boundary elements of all markers. ---*/
  for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
    for(unsigned long l=0; l<nElem_Bound[iMarker]; ++l) {

      /* Get the number of donor elements for the wall function treatment
         and the pointer to the array which stores this info. */
      const unsigned short nDonors = bound[iMarker][l]->GetNDonorsWallFunctions();
      const unsigned long  *donors = bound[iMarker][l]->GetDonorsWallFunctions();

      /*--- Loop over the number of donors for this boundary element. ---*/
      for(unsigned short i=0; i<nDonors; ++i) {

        /*--- Check if the donor element is an external element. ---*/
        unordered_map<unsigned long,unsigned long>::const_iterator UMI;
        UMI = Global_to_Local_Elem.find(donors[i]);

        if(UMI == Global_to_Local_Elem.end()) {

          /*--- Check if element is not already present in
                mapExternalElemIDToTimeLevel. ---*/
          map<unsigned long,CUnsignedShort2T>::iterator MI;
          MI = mapExternalElemIDToTimeLevel.find(donors[i]);
          if(MI == mapExternalElemIDToTimeLevel.end()) {

            /*--- Element not present in external. Add it. ---*/
            mapExternalElemIDToTimeLevel[donors[i]] = CUnsignedShort2T(0,0);
          }

          /*--- The reverse connection may not be present either. Store the global
                ID of this element in the send buffers for the additional
                externals. ---*/
          const unsigned long rankDonor = elemPartitioner.GetRankContainingIndex(donors[i]);

          sendBufAddExternals[rankDonor].push_back(bound[iMarker][l]->GetDomainElement());
          recvFromRank[rankDonor] = 1;
        }
      }
    }
  }

#ifdef HAVE_MPI

  /*--- Determine the number of messages this rank will receive with additional
        externals to be stored. ---*/
  int nRankRecv;
  vector<int> sizeSend(size, 1);
  SU2_MPI::Reduce_scatter(recvFromRank.data(), &nRankRecv, sizeSend.data(),
                          MPI_INT, MPI_SUM, SU2_MPI::GetComm());

  /*--- Determine the number of messages this rank will send. ---*/
  int nRankSend = 0;
  for(int i=0; i<size; ++i)
    nRankSend += recvFromRank[i];

  /*--- Send the data using non-blocking sends to avoid deadlock. ---*/
  vector<SU2_MPI::Request> sendReqs(nRankSend);
  nRankSend = 0;
  for(int i=0; i<size; ++i) {
    if( recvFromRank[i] ) {
      sort(sendBufAddExternals[i].begin(), sendBufAddExternals[i].end());
      vector<unsigned long>::iterator lastElem = unique(sendBufAddExternals[i].begin(),
                                                        sendBufAddExternals[i].end());
      sendBufAddExternals[i].erase(lastElem, sendBufAddExternals[i].end());

      SU2_MPI::Isend(sendBufAddExternals[i].data(), sendBufAddExternals[i].size(),
                     MPI_UNSIGNED_LONG, i, i, SU2_MPI::GetComm(), &sendReqs[nRankSend++]);
    }
  }

  /*--- Loop over the number of ranks from which this rank will receive data
        to be stored in mapExternalElemIDToTimeLevel. ---*/
  for(int i=0; i<nRankRecv; ++i) {

    /*--- Block until a message arrives and determine the source and size
          of the message. Allocate the memory for a receive buffer. ---*/
    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    int source = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);
    vector<unsigned long> recvBuf(sizeMess);

    SU2_MPI::Recv(recvBuf.data(), sizeMess, MPI_UNSIGNED_LONG,
                  source, rank, SU2_MPI::GetComm(), &status);

    /*--- Loop over the entries of recvBuf and add them to
          mapExternalElemIDToTimeLevel, if not present already. ---*/
    for(int j=0; j<sizeMess; ++j) {
      map<unsigned long,CUnsignedShort2T>::iterator MI;
      MI = mapExternalElemIDToTimeLevel.find(recvBuf[j]);
      if(MI == mapExternalElemIDToTimeLevel.end())
        mapExternalElemIDToTimeLevel[recvBuf[j]] = CUnsignedShort2T(0,0);
    }
  }

  /*--- Complete the non-blocking sends. Synchronize the processors afterwards,
        because wild cards have been used in the communication. ---*/
  SU2_MPI::Waitall(nRankSend, sendReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#endif

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Initialize the time level of the owned elements.           ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Get the number of time levels used and check whether or not
        time accurate local time stepping is used. ---*/
  const unsigned short nTimeLevels = config->GetnLevels_TimeAccurateLTS();

  if(nTimeLevels == 1) {

    /*--- No time accurate local time stepping. Set the time level of
          all elements to zero. ---*/
    for(unsigned long i=0; i<nElem; ++i)
      elem[i]->SetTimeLevel(0);
  }
  else {

    /*--- Time accurate local time stepping is used. The estimate of the time
          step is based on free stream values at the moment, but this is easy
          to change, if needed. ---*/
    const su2double Gamma   = config->GetGamma();
    const su2double Prandtl = config->GetPrandtl_Lam();

    const su2double Density   = config->GetDensity_FreeStreamND();
    const su2double *Vel      = config->GetVelocity_FreeStreamND();
    const su2double Viscosity = config->GetViscosity_FreeStreamND();

    su2double VelMag = 0.0;
    for(unsigned short iDim=0; iDim<nDim; ++iDim)
      VelMag += Vel[iDim]*Vel[iDim];
    VelMag = sqrt(VelMag);

    const su2double SoundSpeed = VelMag/config->GetMach();

    /*--- In the estimate of the time step the spectral radius of the inviscid terms is
          needed. As the current estimate is based on the free stream, the value of this
          spectral radius can be computed beforehand. Note that this is a rather
          conservative estimate. ---*/
    su2double charVel2 = 0.0;
    for(unsigned short iDim=0; iDim<nDim; ++iDim) {
      const su2double rad = fabs(Vel[iDim]) + SoundSpeed;
      charVel2 += rad*rad;
    }

    const su2double charVel = sqrt(charVel2);

    /*--- Also the viscous contribution to the time step is constant. Compute it. ---*/
    const su2double factHeatFlux  =  Gamma/Prandtl;
    const su2double lambdaOverMu  = -TWO3;
    const su2double radVisc       =  max(max(1.0, 2.0+lambdaOverMu),factHeatFlux)
                                  *  Viscosity/Density;

    /*--- Allocate the memory for time step estimate of the local elements and
          determine the values. Also keep track of the minimum value. ---*/
    su2double minDeltaT = 1.e25;
    vector<su2double> timeStepElements(nElem);

    for(unsigned long i=0; i<nElem; ++i) {

      unsigned short nPoly = elem[i]->GetNPolySol();
      if(nPoly == 0) nPoly = 1;
      const su2double lenScaleInv = nPoly/elem[i]->GetLengthScale();
      const su2double lenScale    = 1.0/lenScaleInv;

      timeStepElements[i] = lenScale/(charVel + lenScaleInv*radVisc);
      minDeltaT           = min(minDeltaT, timeStepElements[i]);
    }

    /*--- Determine the minimum value of all elements in the grid.
          Only needed for a parallel implementation. ---*/
#ifdef HAVE_MPI
    su2double locVal = minDeltaT;
    SU2_MPI::Allreduce(&locVal, &minDeltaT, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
#endif

    /*--- Initial estimate of the time level of the owned elements. ---*/
    for(unsigned long i=0; i<nElem; ++i) {
      unsigned short timeLevel;
      su2double deltaT = minDeltaT;
      for(timeLevel=0; timeLevel<(nTimeLevels-1); ++timeLevel) {
        deltaT *= 2;
        if(timeStepElements[i] < deltaT) break;
      }

      elem[i]->SetTimeLevel(timeLevel);
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Set up the variables to carry out the MPI communication    ---*/
  /*---         of the external element data.                              ---*/
  /*--------------------------------------------------------------------------*/

  map<unsigned long,CUnsignedShort2T>::iterator MI;
#ifdef HAVE_MPI

  /*--- Determine the ranks from which I receive element data during
        the actual exchange. ---*/
  recvFromRank.assign(size, 0);

  for(MI =mapExternalElemIDToTimeLevel.begin();
      MI!=mapExternalElemIDToTimeLevel.end(); ++MI) {

    /*--- Determine the rank where this external is stored. ---*/
    const unsigned long rankElem = elemPartitioner.GetRankContainingIndex(MI->first);

    /*--- Set the corresponding index of recvFromRank to 1. ---*/
    recvFromRank[rankElem] = 1;
  }

  /*--- Store the ranks from which I receive data in a map. ---*/
  map<int,int> mapRankToIndRecv;
  for(int i=0; i<size; ++i) {
    if( recvFromRank[i] ) {
      int ind = mapRankToIndRecv.size();
      mapRankToIndRecv[i] = ind;
    }
  }

  /*--- Determine the number of ranks from which I will receive data and to
        which I will send data. ---*/
  nRankRecv = mapRankToIndRecv.size();
  SU2_MPI::Reduce_scatter(recvFromRank.data(), &nRankSend, sizeSend.data(),
                          MPI_INT, MPI_SUM, SU2_MPI::GetComm());

  /*--- Create the vector of vectors of the global element ID's that
        will be received from other ranks. ---*/
  vector<vector<unsigned long> > recvElem(nRankRecv, vector<unsigned long>(0));

  for(MI =mapExternalElemIDToTimeLevel.begin();
      MI!=mapExternalElemIDToTimeLevel.end(); ++MI) {

    const unsigned long elemID   = MI->first;
    const unsigned long rankElem = elemPartitioner.GetRankContainingIndex(elemID);

    map<int,int>::const_iterator MRI = mapRankToIndRecv.find(rankElem);
    recvElem[MRI->second].push_back(elemID);
  }

  /*--- Loop over the ranks from which I receive data during the actual
        exchange and send over the global element ID's. In order to avoid
        unnecessary communication, multiple entries are filtered out. ---*/
  sendReqs.resize(nRankRecv);
  map<int,int>::const_iterator MRI = mapRankToIndRecv.begin();

  for(int i=0; i<nRankRecv; ++i, ++MRI) {

    sort(recvElem[i].begin(), recvElem[i].end());
    vector<unsigned long>::iterator lastElem = unique(recvElem[i].begin(),
                                                      recvElem[i].end());
    recvElem[i].erase(lastElem, recvElem[i].end());

    SU2_MPI::Isend(recvElem[i].data(), recvElem[i].size(), MPI_UNSIGNED_LONG,
                   MRI->first, MRI->first, SU2_MPI::GetComm(), &sendReqs[i]);
  }

  /*--- Receive the messages in arbitrary sequence and store the requested
        element ID's, which are converted to local ID's. Furthermore, store
        the processors from which the requested element ID's came from. ---*/
  vector<vector<unsigned long> > sendElem(nRankSend, vector<unsigned long>(0));
  vector<int> sendRank(nRankSend);

  for(int i=0; i<nRankSend; ++i) {

    SU2_MPI::Status status;
    SU2_MPI::Probe(MPI_ANY_SOURCE, rank, SU2_MPI::GetComm(), &status);
    sendRank[i] = status.MPI_SOURCE;

    int sizeMess;
    SU2_MPI::Get_count(&status, MPI_UNSIGNED_LONG, &sizeMess);
    sendElem[i].resize(sizeMess);

    SU2_MPI::Recv(sendElem[i].data(), sizeMess, MPI_UNSIGNED_LONG,
                  sendRank[i], rank, SU2_MPI::GetComm(), &status);

    for(int j=0; j<sizeMess; ++j)
      sendElem[i][j] -= elemPartitioner.GetFirstIndexOnRank(rank);
  }

  /*--- Complete the non-blocking sends. Synchronize the processors afterwards,
        because wild cards have been used in the communication. ---*/
  SU2_MPI::Waitall(nRankRecv, sendReqs.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Barrier(SU2_MPI::GetComm());

#endif

  /*--------------------------------------------------------------------------*/
  /*--- Step 4: Communicate the data of the externals. This data is the    ---*/
  /*---         time level and the number of DOFs of the element.          ---*/
  /*--------------------------------------------------------------------------*/

#ifdef HAVE_MPI

  /*--- Define the send buffers and resize the vector of send requests. ---*/
  vector<vector<unsigned short> > sendBuf(nRankSend, vector<unsigned short>(0));
  sendReqs.resize(nRankSend);

  /*--- Define the return buffers and the vector of return requests. ---*/
  vector<vector<unsigned short> > returnBuf(nRankRecv, vector<unsigned short>(0));
  vector<SU2_MPI::Request> returnReqs(nRankRecv);

  /*--- Copy the information of the time level and number of DOFs into the
        send buffers and send the data using non-blocking sends. ---*/
  for(int i=0; i<nRankSend; ++i) {

    sendBuf[i].resize(2*sendElem[i].size());
    for(unsigned long j=0; j<sendElem[i].size(); ++j) {
      sendBuf[i][2*j]   = elem[sendElem[i][j]]->GetTimeLevel();
      sendBuf[i][2*j+1] = elem[sendElem[i][j]]->GetNDOFsSol();
    }

    SU2_MPI::Isend(sendBuf[i].data(), sendBuf[i].size(), MPI_UNSIGNED_SHORT,
                   sendRank[i], sendRank[i], SU2_MPI::GetComm(), &sendReqs[i]);
  }

  /*--- Receive the data for the externals. As this data is needed immediately,
        blocking communication is used. The time level and the number of DOFs
        of the externals is stored in the second entry of the
        mapExternalElemIDToTimeLevel, which is set accordingly. ---*/
  MRI = mapRankToIndRecv.begin();
  for(int i=0; i<nRankRecv; ++i, ++MRI) {

    returnBuf[i].resize(2*recvElem[i].size());
    SU2_MPI::Status status;
    SU2_MPI::Recv(returnBuf[i].data(), returnBuf[i].size(), MPI_UNSIGNED_SHORT,
                  MRI->first, rank, SU2_MPI::GetComm(), &status);

    for(unsigned long j=0; j<recvElem[i].size(); ++j) {
      MI = mapExternalElemIDToTimeLevel.find(recvElem[i][j]);
      if(MI == mapExternalElemIDToTimeLevel.end())
        SU2_MPI::Error("Entry not found in mapExternalElemIDToTimeLevel", CURRENT_FUNCTION);
      MI->second.short0 = returnBuf[i][2*j];
      MI->second.short1 = returnBuf[i][2*j+1];
    }
  }

  /*--- Complete the nonblocking sends. ---*/
  SU2_MPI::Waitall(nRankSend, sendReqs.data(), MPI_STATUSES_IGNORE);

#endif

  /*--------------------------------------------------------------------------*/
  /*--- Step 5: Only when time accurate local time stepping is employed.   ---*/
  /*---         Iterative algorithm to make sure that neighboring elements ---*/
  /*---         differ at most one time level. This is no fundamental      ---*/
  /*---         issue, but it makes the parallel implementation easier,    ---*/
  /*---         while the penalty in efficiency should be small for most   ---*/
  /*---         cases. Furthermore, also for practical reasons, make sure  ---*/
  /*---         that the donor elements for the wall function treatment    ---*/
  /*---         have the same time level as the element to which the       ---*/
  /*---         boundary face belongs.                                     ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Test for time accurate local time stepping. ---*/
  if(nTimeLevels > 1) {

    /*--- Infinite loop for the iterative algorithm. ---*/
    for(;;) {

      /*--- Variable to indicate whether the local situation has changed. ---*/
      unsigned short localSituationChanged = 0;

      /*--- Loop over the boundary markers and make sure that the time level
            of the element adjacent to the boundary and the donor element for
            the wall function treatment is the same. This is not much of a
            restriction, because in the vast majority of cases, the donor
            element is this adjacent element itself. Requiring it to be the
            same makes the implementation of the wall functions easier. ---*/
      for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
        for(unsigned long l=0; l<nElem_Bound[iMarker]; ++l) {

          /*--- Determine the ID of the adjacent element. ---*/
          const unsigned long elemID = bound[iMarker][l]->GetDomainElement()
                                     - elemPartitioner.GetFirstIndexOnRank(rank);

          /*--- Get the number of donor elements for the wall function treatment
                and the pointer to the array which stores this info. ---*/
          const unsigned short nDonors = bound[iMarker][l]->GetNDonorsWallFunctions();
          const unsigned long  *donors = bound[iMarker][l]->GetDonorsWallFunctions();

          /* Loop over the number of donors and check the time levels. */
          for(unsigned short i=0; i<nDonors; ++i) {

            /*--- Determine the status of the donor element. ---*/
            if(elemPartitioner.GetRankContainingIndex(donors[i]) == static_cast<unsigned long>(rank)) {

              /*--- Donor is stored locally. Determine its local ID and
                    get the time levels of both elements. ---*/
              const unsigned long  donorID    = donors[i]
                                              - elemPartitioner.GetFirstIndexOnRank(rank);
              const unsigned short timeLevelB = elem[elemID]->GetTimeLevel();
              const unsigned short timeLevelD = elem[donorID]->GetTimeLevel();
              const unsigned short timeLevel  = min(timeLevelB, timeLevelD);

              /*--- If the time level of either element is larger than timeLevel,
                    adapt the time levels and indicate that the local situation
                    has changed. ---*/
              if(timeLevelB > timeLevel) {
                elem[elemID]->SetTimeLevel(timeLevel);
                localSituationChanged = 1;
              }

              if(timeLevelD > timeLevel) {
                elem[donorID]->SetTimeLevel(timeLevel);
                localSituationChanged = 1;
              }
            }
            else {

              /*--- The donor element is stored on a different processor.
                    Retrieve its time level from mapExternalElemIDToTimeLevel
                    and determine the minimum time level. ---*/
              const unsigned short timeLevelB = elem[elemID]->GetTimeLevel();
              MI = mapExternalElemIDToTimeLevel.find(donors[i]);
              if(MI == mapExternalElemIDToTimeLevel.end())
                SU2_MPI::Error("Entry not found in mapExternalElemIDToTimeLevel",
                               CURRENT_FUNCTION);
              const unsigned short timeLevel = min(timeLevelB, MI->second.short0);

              /*--- If the time level of either element is larger than timeLevel,
                    adapt the time levels and indicate that the local situation
                    has changed. ---*/
              if(timeLevelB > timeLevel) {
                elem[elemID]->SetTimeLevel(timeLevel);
                localSituationChanged = 1;
              }

              if(MI->second.short0 > timeLevel) {
                MI->second.short0 = timeLevel;
                localSituationChanged = 1;
              }
            }
          }
        }
      }

      /*--- Loop over the matching faces and update the time levels of the
            adjacent elements, if needed. ---*/
      for(vector<CFaceOfElement>::const_iterator FI =localFaces.begin();
                                                 FI!=localFaces.end(); ++FI) {
        /*--- Safeguard against non-matching faces. ---*/
        if(FI->elemID1 < Global_nElem) {

          /*--- Local element ID of the first element. Per definition this is
                always a locally stored element. Also store its time level. ---*/
          const unsigned long  elemID0    = FI->elemID0
                                          - elemPartitioner.GetFirstIndexOnRank(rank);
          const unsigned short timeLevel0 = elem[elemID0]->GetTimeLevel();

          /*--- Determine the status of the second element. ---*/
          if(elemPartitioner.GetRankContainingIndex(FI->elemID1) == static_cast<unsigned long>(rank)) {

            /*--- Both elements are stored locally. Determine the local
                  element of the second element and determine the minimum
                  time level. ---*/
            const unsigned long  elemID1    = FI->elemID1 - elemPartitioner.GetFirstIndexOnRank(rank);
            const unsigned short timeLevel1 = elem[elemID1]->GetTimeLevel();
            const unsigned short timeLevel  = min(timeLevel0, timeLevel1);

            /*--- If the time level of either element is larger than timeLevel+1,
                  adapt the time levels and indicate that the local situation
                  has changed. ---*/
            if(timeLevel0 > timeLevel+1) {
              elem[elemID0]->SetTimeLevel(timeLevel+1);
              localSituationChanged = 1;
            }

            if(timeLevel1 > timeLevel+1) {
              elem[elemID1]->SetTimeLevel(timeLevel+1);
              localSituationChanged = 1;
            }
          }
          else {

            /*--- The second element is stored on a different processor.
                  Retrieve its time level from mapExternalElemIDToTimeLevel
                  and determine the minimum time level. ---*/
            MI = mapExternalElemIDToTimeLevel.find(FI->elemID1);
            if(MI == mapExternalElemIDToTimeLevel.end())
              SU2_MPI::Error("Entry not found in mapExternalElemIDToTimeLevel",
                             CURRENT_FUNCTION);
            const unsigned short timeLevel = min(timeLevel0, MI->second.short0);

            /*--- If the time level of either element is larger than timeLevel+1,
                  adapt the time levels and indicate that the local situation
                  has changed. ---*/
            if(timeLevel0 > timeLevel+1) {
              elem[elemID0]->SetTimeLevel(timeLevel+1);
              localSituationChanged = 1;
            }

            if(MI->second.short0 > timeLevel+1) {
              MI->second.short0 = timeLevel+1;
              localSituationChanged = 1;
            }
          }
        }
      }

      /* Determine whether or not the global situation changed. If not
         the infinite loop can be terminated. */
      unsigned short globalSituationChanged = localSituationChanged;

#ifdef HAVE_MPI
      SU2_MPI::Allreduce(&localSituationChanged, &globalSituationChanged,
                         1, MPI_UNSIGNED_SHORT, MPI_MAX, SU2_MPI::GetComm());
#endif
      if( !globalSituationChanged ) break;

      /*--- Communicate the information of the externals, if needed. ---*/
#ifdef HAVE_MPI

      /*--- Copy the information of the time level into the send buffers
            and send the data using non-blocking sends. Note that the size
            of sendElem is used and not sendBuf, because the size of sendBuf
            is twice the size, see step 4.  ---*/
      for(int i=0; i<nRankSend; ++i) {

        for(unsigned long j=0; j<sendElem[i].size(); ++j)
          sendBuf[i][j] = elem[sendElem[i][j]]->GetTimeLevel();

        SU2_MPI::Isend(sendBuf[i].data(), sendElem[i].size(), MPI_UNSIGNED_SHORT,
                       sendRank[i], sendRank[i], SU2_MPI::GetComm(), &sendReqs[i]);
      }

      /*--- Receive the data for the externals. As this data is needed
            immediately, blocking communication is used. The time level of the
            externals is stored in mapExternalElemIDToTimeLevel, whose second
            element, which contains the time level, is updated accordingly.
            Note that the minimum of the two levels is taken, such that changes
            for the external are also incorporated. As such this information
            must be returned to the original rank. Also note that the size
            of recvElem is used and not returnBuf, because the latter is twice
            as large, see step 4. ---*/
      MRI = mapRankToIndRecv.begin();
      for(int i=0; i<nRankRecv; ++i, ++MRI) {

        SU2_MPI::Status status;
        SU2_MPI::Recv(returnBuf[i].data(), recvElem[i].size(), MPI_UNSIGNED_SHORT,
                      MRI->first, rank, SU2_MPI::GetComm(), &status);

        for(unsigned long j=0; j<recvElem[i].size(); ++j) {
          MI = mapExternalElemIDToTimeLevel.find(recvElem[i][j]);
          if(MI == mapExternalElemIDToTimeLevel.end())
            SU2_MPI::Error("Entry not found in mapExternalElemIDToTimeLevel", CURRENT_FUNCTION);
          MI->second.short0 = min(returnBuf[i][j], MI->second.short0);
          returnBuf[i][j]   = MI->second.short0;
        }

        SU2_MPI::Isend(returnBuf[i].data(), recvElem[i].size(), MPI_UNSIGNED_SHORT,
                       MRI->first, MRI->first+1, SU2_MPI::GetComm(), &returnReqs[i]);
      }

      /*--- Complete the first round of nonblocking sends, such that the
            original send buffers can be used as receive buffers for the
            second round of messages. ---*/
      SU2_MPI::Waitall(nRankSend, sendReqs.data(), MPI_STATUSES_IGNORE);

      /*--- Loop again over the original sending processors to receive
            the updated time level of elements that may have been updated
            on other ranks. ---*/
      for(int i=0; i<nRankSend; ++i) {

        SU2_MPI::Status status;
        SU2_MPI::Recv(sendBuf[i].data(), sendElem[i].size(), MPI_UNSIGNED_SHORT,
                      sendRank[i], rank+1, SU2_MPI::GetComm(), &status);

        for(unsigned long j=0; j<sendElem[i].size(); ++j)
          elem[sendElem[i][j]]->SetTimeLevel(sendBuf[i][j]);
      }

      /* Complete the second round of nonblocking sends. */
      SU2_MPI::Waitall(nRankRecv, returnReqs.data(), MPI_STATUSES_IGNORE);
#endif
    }

    /*------------------------------------------------------------------------*/
    /*--- Step 6: Print a nice output about the number of elements per     ---*/
    /*---         time level.                                              ---*/
    /*------------------------------------------------------------------------*/

    if(rank == MASTER_NODE)
      cout << endl <<"------- Element distribution for time accurate local time stepping ------" << endl;

    /*--- Determine the local number of elements per time level. ---*/
    vector<unsigned long> nLocalElemPerLevel(nTimeLevels, 0);
    for(unsigned long i=0; i<nElem; ++i)
      ++nLocalElemPerLevel[elem[i]->GetTimeLevel()];

    /*--- Determine the global version of nLocalElemPerLevel. This only needs to
          be known on the master node. ---*/
    vector<unsigned long> nGlobalElemPerLevel = nLocalElemPerLevel;

#ifdef HAVE_MPI
     SU2_MPI::Reduce(nLocalElemPerLevel.data(), nGlobalElemPerLevel.data(),
                     nTimeLevels, MPI_UNSIGNED_LONG, MPI_SUM,
                     MASTER_NODE, SU2_MPI::GetComm());
#endif

    /*--- Write the output. ---*/
    if(rank == MASTER_NODE) {
      for(unsigned short i=0; i<nTimeLevels; ++i) {
        if( nGlobalElemPerLevel[i] )
          cout << "Number of elements time level " << i << ": "
               << nGlobalElemPerLevel[i] << endl;
      }
    }

    /*--- And a newline, such that the output is more visible. ---*/
    if(rank == MASTER_NODE) cout << endl;
  }
}

void CPhysicalGeometry::StoreFaceInfoInLocalElements(const vector<CFaceOfElement> &localMatchingFaces) {

  /*--- Define the linear partitioning of the elements. ---*/
  CLinearPartitioner elemPartitioner(Global_nElem, 0);

  /*--- Loop over the locally stored elements. ---*/
  for(unsigned long k=0; k<nElem; ++k) {

    /*--- Get the number of faces and the corner points of these
          faces of this element. ---*/
    unsigned short nFaces;
    unsigned short nPointsPerFace[6];
    unsigned long  faceConn[6][4];

    elem[k]->GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);

    /*--- Initialize the data structures for the neighbors. ---*/
    elem[k]->InitializeNeighbors(nFaces);

    /*--- Loop over the faces of the element. ---*/
    for(unsigned short i=0; i<nFaces; ++i) {

      /*--- Create an object of the class CFaceOfElement for this face. ---*/
      CFaceOfElement thisFace;
      thisFace.nCornerPoints = nPointsPerFace[i];
      for(unsigned short j=0; j<nPointsPerFace[i]; ++j)
        thisFace.cornerPoints[j] = faceConn[i][j];
      thisFace.elemID0    = k + elemPartitioner.GetFirstIndexOnRank(rank);
      thisFace.nPolySol0  = elem[k]->GetNPolySol();
      thisFace.nDOFsElem0 = elem[k]->GetNDOFsSol();

      thisFace.CreateUniqueNumbering();

      /*--- Search for this face in localMatchingFaces. ---*/
      vector<CFaceOfElement>::const_iterator low;
      low = lower_bound(localMatchingFaces.begin(), localMatchingFaces.end(), thisFace);

      /*--- Store the neighboring element and the ownership of this face.
            Store the periodic index for a periodic face. ---*/
      if(low != localMatchingFaces.end() ) {
        if( !(thisFace < *low) ) {
          if(low->elemID0 == thisFace.elemID0) {
            elem[k]->SetNeighbor_Elements(low->elemID1, i);
            elem[k]->SetOwnerFace(low->elem0IsOwner, i);
          }
          else {
            elem[k]->SetNeighbor_Elements(low->elemID0, i);
            elem[k]->SetOwnerFace(!(low->elem0IsOwner), i);
          }

          if(low->periodicIndex > 0)
            elem[k]->SetPeriodicIndex(low->periodicIndex-1, i);
        }
      }
    }
  }
}
