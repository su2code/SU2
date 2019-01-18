/*!
 * \file fem_interpolation_structure.cpp
 * \brief Functions for interpolation for the FEM solver.
 * \author B. Munguía, J. Mukhopadhaya, E. van der Weide
 * \version 6.1.0 "Falcon"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

#include "../include/fem_interpolation_structure.hpp"

CFEMInterpolationSol::CFEMInterpolationSol(void){}

CFEMInterpolationSol::CFEMInterpolationSol(CConfig**      config,
                                           CGeometry**    input_geometry,
                                           CGeometry**    output_geometry,
                                           CSolver***     input_solution,
                                           CSolver***     output_solution,
                                           unsigned short nZone)
{
  // Load geometry data into interpolation grid class
  CFEMInterpolationGrid InputGrid(config, input_geometry, nZone);
  CFEMInterpolationGrid OutputGrid(config, output_geometry, nZone);
  
}

CFEMInterpolationSol::~CFEMInterpolationSol(void){}

void CFEMInterpolationSol::InterpolateSolution(
                                   const std::vector<std::vector<su2double> > &coorInterpol,
                                   const CFEMInterpolationGrid                *inputGrid,
                                   const CFEMInterpolationSol                 *inputSol,
                                   const CFEMInterpolationGrid                *outputGrid)
{
  // Determine the total number of DOFs for which memory must be allocated.
  const unsigned short nZones = inputGrid->GetNZones();
  const unsigned short nDim   = inputGrid->GetNDim();
  unsigned long nDOFsTot = 0;
  for(unsigned short zone=0; zone<nZones; ++zone)
    nDOFsTot += coorInterpol[zone].size();
  nDOFsTot /= nDim;
  
  // Copy the meta data from the input solution.
  for(int i=0; i<5; ++i)
    mHeaderFile[i] = inputSol->mHeaderFile[i];
  mHeaderFile[2] = nDOFsTot;
  
  for(int i=0; i<8; ++i)
    mMetaData[i] = inputSol->mMetaData[i];
  
  mIterNumber = inputSol->mIterNumber;
  mVarNames   = inputSol->mVarNames;
  
  // Determine the number of variables to be interpolated and allocate the memory
  // for the solution DOFs.
  const unsigned short nVar = inputSol->GetNVar();
  
  mSolDOFs.resize(nDOFsTot);
  for(unsigned long l=0; l<nDOFsTot; ++l)
    mSolDOFs[l].resize(nVar);
  
  // Easier storage of the solution format of the input grid.
  const SolutionFormatT solFormatInput = inputGrid->GetSolutionFormat();
  
  // Initialize the zone offset for the input and output solution to zero.
  unsigned long zoneOffsetInputSol  = 0;
  unsigned long zoneOffsetOutputSol = 0;
  
  // Loop over the number of zones.
  for(unsigned short zone=0; zone<nZones; ++zone)
  {
    // Get the zone for the input and output grid as a constant pointer.
    const CFEMInterpolationGridZone *inputGridZone  = inputGrid->GetGridZone(zone);
    const CFEMInterpolationGridZone *outputGridZone = outputGrid->GetGridZone(zone);
    
    // Apply a correction to the coordinates when curved boundaries
    // are present.
    std::vector<su2double> coorInterpolZone;
    ApplyCurvatureCorrection(zone, nDim, inputGridZone, outputGridZone,
                             coorInterpol[zone], coorInterpolZone);
    
    // Define the vectors of the standard elements and the vector to store the
    // standard element for the volume elements.
    std::vector<CFEMStandardElement> standardElementsGrid;
    std::vector<CFEMStandardElement> standardElementsSol;
    std::vector<unsigned short> indInStandardElements;
    
    // First carry out a volume interpolation. Keep track of the points that
    // do not fall within the grid (typically due to a different discrete
    // representation of the boundary of the domain).
    std::vector<unsigned long> pointsForMinDistance;
    VolumeInterpolationSolution(zone, coorInterpolZone, inputGridZone, inputSol,
                                zoneOffsetInputSol, zoneOffsetOutputSol,
                                solFormatInput, pointsForMinDistance,
                                standardElementsGrid, standardElementsSol,
                                indInStandardElements);
    
    // Carry out a surface interpolation, via a minimum distance search,
    // for the points that could not be interpolated via the regular volume
    // interpolation. Print a warning about this.
    if( pointsForMinDistance.size() )
    {
      std::cout << "Zone " << zone << ": " << pointsForMinDistance.size()
      << " DOFs for which the containment search failed." << std::endl;
      std::cout << "A minimum distance search to the boundary of the "
      << "domain is used for these points. " << std::endl;
      
      SurfaceInterpolationSolution(zone, coorInterpolZone, inputGridZone, inputSol,
                                   zoneOffsetInputSol, zoneOffsetOutputSol,
                                   solFormatInput, pointsForMinDistance,
                                   standardElementsSol, indInStandardElements);
    }
    
    // Update the zone offset for the input and output solution.
    zoneOffsetInputSol  += inputGridZone->GetNSolDOFs(solFormatInput);
    zoneOffsetOutputSol += coorInterpol[zone].size()/nDim;
  }
}

void CFEMInterpolationSol::ApplyCurvatureCorrection(
                                        const unsigned short                     zoneID,
                                        const unsigned short                     nDim,
                                        const CFEMInterpolationGridZone          *inputGridZone,
                                        const CFEMInterpolationGridZone          *outputGridZone,
                                        const std::vector<su2double>             &coorOriginal,
                                        std::vector<su2double>                   &coorCorrected)
{
  // Easier storage of the surface elements and coordinates of the grid zones.
  const std::vector<CFEMInterpolationSurfElem> &inputGridSurfElems      = inputGridZone->mSurfElems;
  const std::vector<std::vector<su2double> > &inputGridCoor = inputGridZone->mCoor;
  
  const std::vector<CFEMInterpolationSurfElem> &outputGridSurfElems      = outputGridZone->mSurfElems;
  const std::vector<std::vector<su2double> > &outputGridCoor = outputGridZone->mCoor;
  
  /*--------------------------------------------------------------------------*/
  /*--- Step 1. Build the surface ADTs for both the input grid and the     ---*/
  /*---         output grid.                                               ---*/
  /*--------------------------------------------------------------------------*/
  
  // Write a message that the surface ADT's are built for this zone.
  std::cout << "Grid zone " << zoneID+1
  << ": Building ADTs of the surface grids ...." << std::flush;
  
  // Define the variables needed for the call to BuildSurfaceADT.
  std::vector<unsigned long>  adjElemID;
  std::vector<unsigned short> faceIDInElement;
  
  std::vector<CFEMStandardBoundaryFace> standardBoundaryFacesGrid;
  std::vector<CFEMStandardBoundaryFace> standardBoundaryFacesSol;
  std::vector<unsigned short> inputGridIndInStandardBoundaryFaces;
  
  // Build the surface ADT for the input grid.
  su2_adtElemClass inputGridADT;
  BuildSurfaceADT(inputGridZone, inputGridADT, standardBoundaryFacesGrid,
                  standardBoundaryFacesSol, inputGridIndInStandardBoundaryFaces,
                  adjElemID, faceIDInElement);
  
  // Build the surface ADT for the output grid.
  std::vector<unsigned short> outputGridIndInStandardBoundaryFaces;
  su2_adtElemClass outputGridADT;
  BuildSurfaceADT(outputGridZone, outputGridADT, standardBoundaryFacesGrid,
                  standardBoundaryFacesSol, outputGridIndInStandardBoundaryFaces,
                  adjElemID, faceIDInElement);
  
  // Write a message that the ADT's were built.
  std::cout << " Done" << std::endl << std::flush;
  
  /*--------------------------------------------------------------------------*/
  /*--- Step 2. Carry out the wall distance searches for both the input    ---*/
  /*---         and the output grid and apply the curvature correction.    ---*/
  /*--------------------------------------------------------------------------*/
  
  // Initialize the corrected coordinates.
  coorCorrected = coorOriginal;
  
  // Determine the write frequency.
  const unsigned long nDOFs = coorOriginal.size()/nDim;
  
  unsigned long writeFreq = nDOFs/5;
  if(writeFreq > 10000) writeFreq = 10000;
  if(writeFreq <   100) writeFreq =   100;
  
  writeFreq = writeFreq/10;
  if(writeFreq == 0) writeFreq = 1;
  writeFreq *= 10;
    
    // Define the vectors used in the tree search. Pre-allocate some memory
    // for efficiency reasons.
    std::vector<su2_BBoxTargetClass> BBoxTargets(200);
    std::vector<unsigned long> frontLeaves(200), frontLeavesNew(200);
    
    // Loop over the DOFs to be corrected.
    for(unsigned long l=0; l<nDOFs; ++l)
    {
      // Write a message about the number of DOFs to be interpolated.
      if( !(l%writeFreq) )
        std::cout << "Grid zone " << zoneID+1 << ": " << l << " out of "
        << nDOFs << " DOFs corrected for surface curvature." << std::endl << std::flush;
      
      // Set a pointer to the coordinates to be searched.
      su2double *coor = coorCorrected.data() + l*nDim;
      
      // Carry out the minimum distance search to the input grid.
      unsigned short subElem;
      unsigned long  parElem;
      int            rank;
      su2double      dist;
      su2double      weightsInterpol[4];
      inputGridADT.DetermineNearestElement(coor, dist, subElem, parElem, rank,
                                           weightsInterpol, BBoxTargets,
                                           frontLeaves, frontLeavesNew);
      
      // Subelement found that minimizes the distance to the given coordinate.
      // However, what is needed is the location in the high order parent element.
      // Determine this.
      su2double parCoor[3], wallCoorInputGrid[3];
      unsigned short ind = inputGridIndInStandardBoundaryFaces[parElem];
      HighOrderMinDistanceSearch(coor, parElem, subElem, weightsInterpol,
                                 &standardBoundaryFacesGrid[ind],
                                 &inputGridSurfElems[parElem],
                                 inputGridCoor, parCoor, wallCoorInputGrid);
      
      // Carry out the minimum distance search to the output grid.
      outputGridADT.DetermineNearestElement(coor, dist, subElem, parElem, rank,
                                            weightsInterpol, BBoxTargets,
                                            frontLeaves, frontLeavesNew);
      
      // Subelement found that minimizes the distance to the given coordinate.
      // However, what is needed is the location in the high order parent element.
      // Determine this.
      su2double wallCoorOutputGrid[3];
      ind = outputGridIndInStandardBoundaryFaces[parElem];
      HighOrderMinDistanceSearch(coor, parElem, subElem, weightsInterpol,
                                 &standardBoundaryFacesGrid[ind],
                                 &outputGridSurfElems[parElem],
                                 outputGridCoor, parCoor, wallCoorOutputGrid);
      
      // Determine the curvature correction, which is the vector from the
      // wall coordinate of the output grid to the wall coordinates on the
      // input grid.
      for(unsigned short iDim=0; iDim<nDim; ++iDim)
        coor[iDim] += wallCoorInputGrid[iDim] - wallCoorOutputGrid[iDim];
    }
  
  // Write a message that the curvature correction is finished.
  std::cout << "Grid zone " << zoneID+1 << ": Curvature correction finished"
  << std::endl << std::flush;
}

void CFEMInterpolationSol::BuildSurfaceADT(
                               const CFEMInterpolationGridZone           *gridZone,
                               su2_adtElemClass                          &surfaceADT,
                               std::vector<CFEMStandardBoundaryFace>     &standardBoundaryFacesGrid,
                               std::vector<CFEMStandardBoundaryFace>     &standardBoundaryFacesSol,
                               std::vector<unsigned short>               &indInStandardBoundaryFaces,
                               std::vector<unsigned long>                &adjElemID,
                               std::vector<unsigned short>               &faceIDInElement)
{
  // Easier storage of the volume elements, surface elements and coordinates of the grid zone.
  const std::vector<CFEMInterpolationVolElem>  &volElems           = gridZone->mVolElems;
  const std::vector<CFEMInterpolationSurfElem> &surfElems          = gridZone->mSurfElems;
  const std::vector<std::vector<su2double> > &coorGrid = gridZone->mCoor;
  const unsigned short nDim    = coorGrid.size();
  const unsigned long  nPoints = coorGrid[0].size();
  
  /*--------------------------------------------------------------------------*/
  /*--- Step 1. Determine the volume elements adjacent to the boundary        */
  /*---         surface elements.                                             */
  /*--------------------------------------------------------------------------*/
  
  // Define the vector, to store the local faces.
  std::vector<CFEMInterpolationFaceOfElem> localFaces;
  
  // Loop over the volume elements to build the vector of faces.
  for(unsigned long l=0; l<volElems.size(); ++l)
  {
    // Determine the corner points of all the faces of this element.
    unsigned short nFaces;
    unsigned short nPointsPerFace[6];
    unsigned long  faceConn[6][4];
    
    volElems[l].GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);
    
    // Loop over the faces, set the appropriate information and store
    // it in localFaces.
    for(unsigned short i=0; i<nFaces; ++i)
    {
      CFEMInterpolationFaceOfElem thisFace;
      thisFace.nCornerPoints = nPointsPerFace[i];
      for(unsigned short j=0; j<nPointsPerFace[i]; ++j)
        thisFace.cornerPoints[j] = faceConn[i][j];
      
      thisFace.elemID = l;
      thisFace.faceID = i;
      
      thisFace.CreateUniqueNumbering();
      localFaces.push_back(thisFace);
    }
  }
  
  // Sort localFaces in increasing order.
  std::sort(localFaces.begin(), localFaces.end());
  
  // Loop over the surface elements to determine the adjacent element and
  // the face ID inside this element.
  adjElemID.resize(surfElems.size());
  faceIDInElement.resize(surfElems.size());
  
  for(unsigned long l=0; l<surfElems.size(); ++l)
  {
    // Determine the corner points of this element.
    unsigned short nPointsPerFace;
    unsigned long  faceConn[4];
    surfElems[l].GetCornerPoints(nPointsPerFace, faceConn);
    
    // Create an object of CFEMInterpolationFaceOfElem to store this info.
    CFEMInterpolationFaceOfElem thisFace;
    thisFace.nCornerPoints = nPointsPerFace;
    for(unsigned short j=0; j<nPointsPerFace; ++j)
      thisFace.cornerPoints[j] = faceConn[j];
    
    thisFace.elemID = 0;
    thisFace.faceID = 0;
    
    thisFace.CreateUniqueNumbering();
    
    // Search for thisFace in localFaces. It must be found.
    if( std::binary_search(localFaces.begin(), localFaces.end(), thisFace) )
    {
      std::vector<CFEMInterpolationFaceOfElem>::const_iterator low;
      low = std::lower_bound(localFaces.begin(), localFaces.end(), thisFace);
      
      // Store the information in adjElemID and faceIDInElement.
      adjElemID[l]       = low->elemID;
      faceIDInElement[l] = low->faceID;
    }
    else
      SU2_MPI::Error("Boundary face not found in volume elements. The grid is not valid.", CURRENT_FUNCTION);
  }
  
  // Release the memory of localFaces again, because it is not needed anymore.
  std::vector<CFEMInterpolationFaceOfElem>().swap(localFaces);
  
  /*--------------------------------------------------------------------------*/
  /*--- Step 2. Build the actual ADT of the surface elements.                 */
  /*--------------------------------------------------------------------------*/
  
  // Initialize an array for the mesh points, which eventually contains the
  // mapping from the local nodes to the number used in the connectivity of the
  // boundary faces. However, in a first pass it is an indicator whether
  // or not a mesh point is on a boundary.
  std::vector<unsigned long> meshToSurface(nPoints, 0);
  
  // Define the vectors, which store the mapping from the subface to the
  // parent face, subface ID within the parent face, the face type and
  // the connectivity of the subfaces.
  std::vector<unsigned long>  parentFace;
  std::vector<unsigned short> subFaceIDInParent;
  std::vector<unsigned short> VTK_TypeFace;
  std::vector<unsigned long>  faceConn;
  
  // Loop over the surface elements to create the connectivity of the subelements.
  indInStandardBoundaryFaces.resize(surfElems.size());
  for(unsigned long l=0; l<surfElems.size(); ++l)
  {
    // Easier storage of the ID of the adjacent element.
    const unsigned long elemID = adjElemID[l];
    
    // Determine the index in the standard boundary faces.
    // If not present yet, create the standard boundary faces.
    unsigned long ind;
    for(ind=0; ind<standardBoundaryFacesGrid.size(); ++ind)
    {
      if(standardBoundaryFacesSol[ind].SameStandardBoundaryFace(surfElems[l].mVTK_TYPE,
                                                                false,
                                                                volElems[elemID].mVTK_TYPE,
                                                                volElems[elemID].mNPolySol,
                                                                false) &&
         standardBoundaryFacesGrid[ind].SameStandardBoundaryFace(surfElems[l].mVTK_TYPE,
                                                                 false,
                                                                 volElems[elemID].mVTK_TYPE,
                                                                 volElems[elemID].mNPolyGrid,
                                                                 false))
        break;
    }
    
    if(ind == standardBoundaryFacesSol.size())
    {
      standardBoundaryFacesSol.push_back(CFEMStandardBoundaryFace(surfElems[l].mVTK_TYPE,
                                                                      volElems[elemID].mVTK_TYPE,
                                                                      volElems[elemID].mNPolySol,
                                                                      false, false) );
      standardBoundaryFacesGrid.push_back(CFEMStandardBoundaryFace(surfElems[l].mVTK_TYPE,
                                                                       volElems[elemID].mVTK_TYPE,
                                                                       volElems[elemID].mNPolyGrid,
                                                                       false, false) );
    }
    
    indInStandardBoundaryFaces[l] = ind;
    
    // Abbreviate the grid DOFs of this element a bit easier.
    const unsigned long *DOFs = surfElems[l].mConnGrid.data();
    
    // Set the flag of the mesh points to true.
    for(unsigned short j=0; j<surfElems[l].mNDOFsGrid; ++j)
      meshToSurface[DOFs[j]] = 1;
    
    // Get the required data from the standard element.
    const unsigned short VTK_Type      = standardBoundaryFacesGrid[ind].GetVTK_Type();
    const unsigned short nSubFaces     = standardBoundaryFacesGrid[ind].GetNSubFaces();
    const unsigned short nDOFsPerFace  = standardBoundaryFacesGrid[ind].GetNDOFsPerSubFace();
    const unsigned short *connSubFaces = standardBoundaryFacesGrid[ind].GetSubFaceConn();
    
    // Loop over the number of subelements and store the required data.
    unsigned short kk = 0;
    for(unsigned short j=0; j<nSubFaces; ++j)
    {
      parentFace.push_back(l);
      subFaceIDInParent.push_back(j);
      VTK_TypeFace.push_back(VTK_Type);
      
      for(unsigned short k=0; k<nDOFsPerFace; ++k, ++kk)
        faceConn.push_back(DOFs[connSubFaces[kk]]);
    }
  }
  
  // Create the coordinates of the points on the surfaces and create the final
  // version of the mapping from all volume points to the surface points.
  std::vector<su2double> surfaceCoor;
  unsigned long nSurfacePoints = 0;
  for(unsigned long i=0; i<nPoints; ++i)
  {
    if( meshToSurface[i] )
    {
      meshToSurface[i] = nSurfacePoints++;
      for(unsigned short k=0; k<nDim; ++k)
        surfaceCoor.push_back(coorGrid[k][i]);
    }
  }
  
  // Change the surface connectivity, such that it corresponds to
  // the entries in surfaceCoor rather than in the volume coordinates.
  for(unsigned long i=0; i<faceConn.size(); ++i)
    faceConn[i] = meshToSurface[faceConn[i]];
  
  // Build the local ADT.
  surfaceADT.CreateADT(nDim, surfaceCoor, faceConn, VTK_TypeFace,
                       subFaceIDInParent, parentFace);
}

void CFEMInterpolationSol::VolumeInterpolationSolution(
                                           const unsigned short                 zoneID,
                                           const std::vector<su2double>         &coorInterpol,
                                           const CFEMInterpolationGridZone      *gridZone,
                                           const CFEMInterpolationSol           *inputSol,
                                           const unsigned long                  zoneOffsetInputSol,
                                           const unsigned long                  zoneOffsetOutputSol,
                                           const SolutionFormatT                solFormatInput,
                                           std::vector<unsigned long>           &pointsSearchFailed,
                                           std::vector<CFEMStandardElement> &standardElementsGrid,
                                           std::vector<CFEMStandardElement> &standardElementsSol,
                                           std::vector<unsigned short>          &indInStandardElements)
{
  /*--------------------------------------------------------------------------*/
  /*--- Step 1. Build the local ADT of the volume elements. Note that the  ---*/
  /*---         ADT is built with the linear subelements. This is done to  ---*/
  /*---         avoid relatively many expensive Newton solves for high     ---*/
  /*---         order elements.                                            ---*/
  /*--------------------------------------------------------------------------*/
  
  // Write a message that the volume ADT is built for this zone.
  std::cout << "Grid zone " << zoneID+1 << ": Building ADT of the volume grid...."
  << std::flush;
  
  // Easier storage of the volume elements and coordinates of the grid zone.
  const std::vector<CFEMInterpolationVolElem> &volElems            = gridZone->mVolElems;
  const std::vector<std::vector<su2double> > &coorGrid = gridZone->mCoor;
  const unsigned short nDim = coorGrid.size();
  
  // Allocate the memory for the vector, which stores the index in the standard
  // elements for the volume elements.
  indInStandardElements.resize(volElems.size());
  
  // Define the vectors, which store the mapping from the subelement to the
  // parent element, subelement ID within the parent element, the element
  // type and the connectivity of the subelements.
  std::vector<unsigned long>  parentElement;
  std::vector<unsigned short> subElementIDInParent;
  std::vector<unsigned short> VTK_TypeElem;
  std::vector<unsigned long>  elemConn;
  
  // Loop over the volume elements to create the connectivity of the subelements.
  for(unsigned long l=0; l<volElems.size(); ++l)
  {
    // Determine the index in the standard elements.
    // If not present yet, create the standard elements.
    unsigned long ind;
    for(ind=0; ind<standardElementsSol.size(); ++ind)
    {
      if(standardElementsSol[ind].SameStandardElement(volElems[l].mVTK_TYPE,
                                                      volElems[l].mNPolySol,
                                                      false) &&
         standardElementsGrid[ind].SameStandardElement(volElems[l].mVTK_TYPE,
                                                       volElems[l].mNPolyGrid,
                                                       false))
        break;
    }
    
    if(ind == standardElementsSol.size())
    {
      standardElementsSol.push_back(CFEMStandardElement(volElems[l].mVTK_TYPE,
                                                            volElems[l].mNPolySol,
                                                            false));
      standardElementsGrid.push_back(CFEMStandardElement(volElems[l].mVTK_TYPE,
                                                             volElems[l].mNPolyGrid,
                                                             false,
                                                             standardElementsSol[ind].GetOrderExact(),
                                                             standardElementsSol[ind].GetRDOFs(),
                                                             standardElementsSol[ind].GetSDOFs(),
                                                             standardElementsSol[ind].GetTDOFs()));
    }
    
    indInStandardElements[l] = ind;
    
    // Determine the necessary data from the corresponding standard element.
    unsigned short VTK_Type[]  = {standardElementsGrid[ind].GetVTK_Type1(),
      standardElementsGrid[ind].GetVTK_Type2()};
    unsigned short nSubElems[] = {0, 0};
    unsigned short nDOFsPerSubElem[] = {0, 0};
    const unsigned short *connSubElems[] = {NULL, NULL};
    
    if(VTK_Type[0] != NONE) {
      nSubElems[0]       = standardElementsGrid[ind].GetNSubElemsType1();
      nDOFsPerSubElem[0] = standardElementsGrid[ind].GetNDOFsPerSubElem(VTK_Type[0]);
      connSubElems[0]    = standardElementsGrid[ind].GetSubConnType1();
    }
    
    if(VTK_Type[1] != NONE) {
      nSubElems[1]       = standardElementsGrid[ind].GetNSubElemsType2();
      nDOFsPerSubElem[1] = standardElementsGrid[ind].GetNDOFsPerSubElem(VTK_Type[1]);
      connSubElems[1]    = standardElementsGrid[ind].GetSubConnType2();
    }
    
    // Abbreviate the grid DOFs of this element a bit easier.
    const unsigned long *DOFs = volElems[l].mConnGrid.data();
    
    // Loop over the number of subelements and store the required data.
    unsigned short jj = 0;
    for(unsigned short i=0; i<2; ++i) {
      unsigned short kk = 0;
      for(unsigned short j=0; j<nSubElems[i]; ++j, ++jj) {
        parentElement.push_back(l);
        subElementIDInParent.push_back(jj);
        VTK_TypeElem.push_back(VTK_Type[i]);
        
        for(unsigned short k=0; k<nDOFsPerSubElem[i]; ++k, ++kk)
          elemConn.push_back(DOFs[connSubElems[i][kk]]);
      }
    }
  }
  
  // Copy the coordinates in a vector that can be used by the ADT.
  std::vector<su2double> volCoor;
  volCoor.reserve(nDim*coorGrid[0].size());
  for(unsigned long l=0; l<coorGrid[0].size(); ++l)
  {
    for(unsigned short k=0; k<nDim; ++k)
      volCoor.push_back(coorGrid[k][l]);
  }
  
  // Build the local ADT.
  su2_adtElemClass volumeADT(nDim, volCoor, elemConn, VTK_TypeElem,
                             subElementIDInParent, parentElement);
  
  // Release the memory of the vectors used to build the ADT. To make sure
  // that all the memory is deleted, the swap function is used.
  std::vector<unsigned short>().swap(subElementIDInParent);
  std::vector<unsigned short>().swap(VTK_TypeElem);
  std::vector<unsigned long>().swap(parentElement);
  std::vector<unsigned long>().swap(elemConn);
  std::vector<su2double>().swap(volCoor);
  
  // Write a message that the ADT was built.
  std::cout << " Done" << std::endl << std::flush;
  
  /*--------------------------------------------------------------------------*/
  /*--- Step 2. Search for donor elements for the given coordinates.       ---*/
  /*--------------------------------------------------------------------------*/
  
  // Determine the write frequency.
  const unsigned long nDOFsInterpol = coorInterpol.size()/nDim;
  
  unsigned long writeFreq = nDOFsInterpol/5;
  if(writeFreq > 10000) writeFreq = 10000;
  if(writeFreq <   100) writeFreq =   100;
  
  writeFreq = writeFreq/10;
  if(writeFreq == 0) writeFreq = 1;
  writeFreq *= 10;
    
    // Define the local vector to store the failed searches for each thread.
    std::vector<long> localPointsFailed;
    
    // Define the vectors used in the tree search. Pre-allocate some memory
    // for efficiency reasons.
    std::vector<unsigned long> frontLeaves(200), frontLeavesNew(200);
    
    // Loop over the DOFs to be interpolated.
    for(unsigned long l=0; l<nDOFsInterpol; ++l)
    {
      // Write a message about the number of DOFs to be interpolated.
      if( !(l%writeFreq) )
        std::cout << "Grid zone " << zoneID+1 << ": " << l << " out of "
        << nDOFsInterpol << " DOFs interpolated." << std::endl << std::flush;
      
      // Set a pointer to the coordinates to be searched.
      const su2double *coor = coorInterpol.data() + l*nDim;
      
      // Carry out the containment search and check if it was successful.
      unsigned short subElem;
      unsigned long  parElem;
      int            rank;
      su2double      parCoor[3], weightsInterpol[8];
      if( volumeADT.DetermineContainingElement(coor, subElem, parElem, rank,
                                               parCoor, weightsInterpol,
                                               frontLeaves, frontLeavesNew) )
      {
        // Subelement found that contains the exchange location. However,
        // what is needed is the location in the high order parent element.
        // Determine this.
        HighOrderContainmentSearch(coor, parElem, subElem, weightsInterpol,
                                   &standardElementsGrid[indInStandardElements[parElem]],
                                   &volElems[parElem], coorGrid, parCoor);
        
        // Carry out the actual interpolation.
        const unsigned short ind = indInStandardElements[parElem];
        const unsigned long  ll  = l + zoneOffsetOutputSol;
        
        SolInterpolate(&standardElementsSol[ind], inputSol, zoneOffsetInputSol,
                       &volElems[parElem], solFormatInput, parCoor, mSolDOFs[ll]);
      }
      else
      {
        // Containment search was not successful. Store the ID in localPointsFailed.
        localPointsFailed.push_back(l);
      }
    }
    
    // Store the local failed points in pointsSearchFailed.
    pointsSearchFailed.insert(pointsSearchFailed.end(),
                              localPointsFailed.begin(), localPointsFailed.end());
  
  // Write a message that the volume search is finished.
  std::cout << "Grid zone " << zoneID+1 << ": Volume search finished"
  << std::endl << std::flush;
}

void CFEMInterpolationSol::HighOrderContainmentSearch(
                                          const su2double                            *coor,
                                          const unsigned long                        parElem,
                                          const unsigned short                       subElem,
                                          const su2double                            *weightsSubElem,
                                          CFEMStandardElement                        *standardElementGrid,
                                          const CFEMInterpolationVolElem             *volElem,
                                          const std::vector<std::vector<su2double> > &coorGrid,
                                          su2double                                  *parCoor)
{
  // Easier storage of the number of dimensions.
  const unsigned short nDim = coorGrid.size();
  
  // Definition of the maximum number of iterations in the Newton solver
  // and the tolerance level. */
  const unsigned short maxIt = 50;
  const su2double tolNewton  = 1.e-10;
  
  /*--------------------------------------------------------------------------*/
  /* Step 1: Create an initial guess for the parametric coordinates from      */
  /*         interpolation in the linear sub-element of the parent element.   */
  /*--------------------------------------------------------------------------*/
  
  // Define the variables to store the number of DOFs and the connectivity
  // of the sub element in which the given coordinate resides.
  unsigned short nDOFsPerSubElem = 0;
  const unsigned short *connSubElems;
  
  // Check if the sub element is of the first sub-element type.
  const unsigned short nSubElemType1 = standardElementGrid->GetNSubElemsType1();
  if(subElem < nSubElemType1) {
    
    // Determine the element type and set nDOFsPerSubElem.
    switch( standardElementGrid->GetVTK_Type1() ) {
      case TRIANGLE:      nDOFsPerSubElem = 3; break;
      case QUADRILATERAL: nDOFsPerSubElem = 4; break;
      case TETRAHEDRON:   nDOFsPerSubElem = 4; break;
      case PYRAMID:       nDOFsPerSubElem = 5; break;
      case PRISM:         nDOFsPerSubElem = 6; break;
      case HEXAHEDRON:    nDOFsPerSubElem = 8; break;
      default: break; // Just to avoid a compiler warning.
    }
    
    // Set the connectivity for the correct subelement.
    connSubElems = standardElementGrid->GetSubConnType1()
    + subElem*nDOFsPerSubElem;
  }
  else {
    
    // The sub-element is of the second sub-element type. Determine the
    // element type and set nDOFsPerSubElem.
    switch( standardElementGrid->GetVTK_Type2() ) {
      case TRIANGLE:      nDOFsPerSubElem = 3; break;
      case QUADRILATERAL: nDOFsPerSubElem = 4; break;
      case TETRAHEDRON:   nDOFsPerSubElem = 4; break;
      case PYRAMID:       nDOFsPerSubElem = 5; break;
      case PRISM:         nDOFsPerSubElem = 6; break;
      case HEXAHEDRON:    nDOFsPerSubElem = 8; break;
      default: break; // Just to avoid a compiler warning.
    }
    
    // Set the connectivity for the correct subelement.
    connSubElems = standardElementGrid->GetSubConnType2()
    + (subElem-nSubElemType1)*nDOFsPerSubElem;
  }
  
  // Get the parametric coordinates of the DOFs from the standard element.
  const std::vector<su2double> *locDOFs[] = {standardElementGrid->GetRDOFs(),
    standardElementGrid->GetSDOFs(),
    standardElementGrid->GetTDOFs()};
  
  /* Create the initial guess of the parametric coordinates by interpolation
   in the sub-element. */
  for(unsigned short iDim=0; iDim<nDim; ++iDim) {
    parCoor[iDim] = 0.0;
    const su2double *coorDOFs = locDOFs[iDim]->data();
    for(unsigned short i=0; i<nDOFsPerSubElem; ++i)
      parCoor[iDim] += weightsSubElem[i]*coorDOFs[connSubElems[i]];
  }
  
  /*--------------------------------------------------------------------------*/
  /* Step 2: The Newton algorithm to compute the parametric coordinates in    */
  /*         the high order element.                                          */
  /*--------------------------------------------------------------------------*/
  
  // Determine the number of DOFs of the standard element and allocate the
  // memory to store the basis functions and its derivatives.
  const unsigned short nDOFs = standardElementGrid->GetNDOFs();
  
  std::vector<su2double> lagBasis(nDOFs);
  std::vector<std::vector<su2double> > dLagBasis(nDim, std::vector<su2double>(nDOFs));
  
  // Abbreviate the grid DOFs of this element a bit easier.
  const unsigned long *DOFs = volElem->mConnGrid.data();
  
  // Loop over the maximum number of iterations.
  unsigned short itCount;
  for(itCount=0; itCount<maxIt; ++itCount) {
    
    // Compute the Lagrangian basis functions and its derivatives in
    // the current parametric coordinate.
    standardElementGrid->BasisFunctionsAndDerivativesInPoint(parCoor, lagBasis,
                                                             dLagBasis);
    
    // Make a distinction between 2D and 3D in order to have the most
    // efficient code.
    bool converged = false;
    switch( nDim ) {
      case 2: {
        // Two dimensional computation. Compute the values of the function
        // and minus the Jacobian matrix.
        su2double f0 = coor[0], f1 = coor[1];
        su2double a00 = 0.0, a01 = 0.0, a10 = 0.0, a11 = 0.0;
        for(unsigned short i=0; i<nDOFs; ++i) {
          const su2double x = coorGrid[0][DOFs[i]];
          const su2double y = coorGrid[1][DOFs[i]];
          
          f0 -= x*lagBasis[i]; f1 -= y*lagBasis[i];
          
          a00 += x*dLagBasis[0][i]; a01 += x*dLagBasis[1][i];
          a10 += y*dLagBasis[0][i]; a11 += y*dLagBasis[1][i];
        }
        
        // Compute the updates of the parametric values. As minus the
        // Jacobian is computed, the updates should be added to parCoor.
        const su2double detInv = 1.0/(a00*a11 - a01*a10);
        const su2double dr = detInv*(f0*a11 - f1*a01);
        const su2double ds = detInv*(f1*a00 - f0*a10);
        
        parCoor[0] += dr;
        parCoor[1] += ds;
        
        // Check for convergence.
        if(fabs(dr) <= tolNewton && fabs(ds) <= tolNewton) converged = true;
        break;
      }
        
      case 3: {
        // Three dimensional computation. Compute the values of the function
        // and minus the Jacobian matrix.
        su2double f0 = coor[0], f1 = coor[1], f2 = coor[2];
        su2double a00 = 0.0, a01 = 0.0, a02 = 0.0;
        su2double a10 = 0.0, a11 = 0.0, a12 = 0.0;
        su2double a20 = 0.0, a21 = 0.0, a22 = 0.0;
        
        for(unsigned short i=0; i<nDOFs; ++i) {
          const su2double x = coorGrid[0][DOFs[i]];
          const su2double y = coorGrid[1][DOFs[i]];
          const su2double z = coorGrid[2][DOFs[i]];
          
          f0 -= x*lagBasis[i]; f1 -= y*lagBasis[i]; f2 -= z*lagBasis[i];
          
          a00 += x*dLagBasis[0][i]; a01 += x*dLagBasis[1][i]; a02 += x*dLagBasis[2][i];
          a10 += y*dLagBasis[0][i]; a11 += y*dLagBasis[1][i]; a12 += y*dLagBasis[2][i];
          a20 += z*dLagBasis[0][i]; a21 += z*dLagBasis[1][i]; a22 += z*dLagBasis[2][i];
        }
        
        // Compute the updates of the parametric values. As minus the
        // Jacobian is computed, the updates should be added to parCoor.
        const su2double detInv = 1.0/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22
                                      +      a01*a12*a20 + a02*a10*a21 - a02*a11*a20);
        const su2double dr =  detInv*(a01*a12*f2 - a01*a22*f1 - a02*a11*f2
                                      +          a02*a21*f1 + a11*a22*f0 - a12*a21*f0);
        const su2double ds = -detInv*(a00*a12*f2 - a00*a22*f1 - a02*a10*f2
                                      +          a02*a20*f1 + a10*a22*f0 - a12*a20*f0);
        const su2double dt =  detInv*(a00*a11*f2 - a00*a21*f1 - a01*a10*f2
                                      +          a01*a20*f1 + a10*a21*f0 - a11*a20*f0);
        parCoor[0] += dr;
        parCoor[1] += ds;
        parCoor[2] += dt;
        
        // Check for convergence.
        if(fabs(dr) <= tolNewton && fabs(ds) <= tolNewton && fabs(dt) <= tolNewton)
          converged = true;
        break;
      }
    }
    
    // Break the loop if the Newton algorithm converged.
    if( converged ) break;
  }
  
  // Terminate if the Newton algorithm did not converge.
  if(itCount == maxIt)
    SU2_MPI::Error("Newton did not converge.", CURRENT_FUNCTION);
}

void CFEMInterpolationSol::SurfaceInterpolationSolution(
                                            const unsigned short                 zoneID,
                                            const std::vector<su2double>         &coorInterpol,
                                            const CFEMInterpolationGridZone      *gridZone,
                                            const CFEMInterpolationSol           *inputSol,
                                            const unsigned long                  zoneOffsetInputSol,
                                            const unsigned long                  zoneOffsetOutputSol,
                                            const SolutionFormatT                solFormatInput,
                                            const std::vector<unsigned long>     &pointsMinDistSearch,
                                            std::vector<CFEMStandardElement>     &standardElementsSol,
                                            const std::vector<unsigned short>    &indInStandardElements)
{
  // Easier storage of the volume elements, the surface elements
  // and coordinates of the grid zone.
  const std::vector<CFEMInterpolationVolElem>  &volElems           = gridZone->mVolElems;
  const std::vector<CFEMInterpolationSurfElem> &surfElems          = gridZone->mSurfElems;
  const std::vector<std::vector<su2double> > &coorGrid = gridZone->mCoor;
  const unsigned short nDim = coorGrid.size();
  
  /*--------------------------------------------------------------------------*/
  /*--- Step 1. Build the local ADT of the surface elements. Note that the ---*/
  /*---         ADT is built with the linear subelements. This is done to  ---*/
  /*---         avoid relatively many expensive Newton solves for high     ---*/
  /*---         order elements.                                            ---*/
  /*--------------------------------------------------------------------------*/
  
  // Write a message that the surface ADT is built for this zone.
  std::cout << "Grid zone " << zoneID+1 << ": Building ADT of the surface grid...."
  << std::flush;
  
  // Define the vectors to store the adjacent element and the face ID
  // inside the element.
  std::vector<unsigned long>  adjElemID;
  std::vector<unsigned short> faceIDInElement;
  
  // Define the vectors of the standard boundary faces and the vector to store
  // the standard boundary face for the surface elements.
  std::vector<CFEMStandardBoundaryFace> standardBoundaryFacesGrid;
  std::vector<CFEMStandardBoundaryFace> standardBoundaryFacesSol;
  std::vector<unsigned short> indInStandardBoundaryFaces;
  
  // Build the surface ADT.
  su2_adtElemClass surfaceADT;
  BuildSurfaceADT(gridZone, surfaceADT, standardBoundaryFacesGrid,
                  standardBoundaryFacesSol, indInStandardBoundaryFaces,
                  adjElemID, faceIDInElement);
  
  // Write a message that the ADT was built.
  std::cout << " Done" << std::endl << std::flush;
  
  /*--------------------------------------------------------------------------*/
  /*--- Step 2. Search for donor elements for the given coordinates.       ---*/
  /*--------------------------------------------------------------------------*/
  
  // Determine the write frequency.
  unsigned long writeFreq = pointsMinDistSearch.size()/5;
  if(writeFreq > 10000) writeFreq = 10000;
  if(writeFreq <   100) writeFreq =   100;
  
  writeFreq = writeFreq/10;
  if(writeFreq == 0) writeFreq = 1;
  writeFreq *= 10;
    
    // Define the vectors used in the tree search. Pre-allocate some memory
    // for efficiency reasons.
    std::vector<su2_BBoxTargetClass> BBoxTargets(200);
    std::vector<unsigned long> frontLeaves(200), frontLeavesNew(200);
    
    // Loop over the points for which a minimum distance search must be carried out.
    for(unsigned long l=0; l<pointsMinDistSearch.size(); ++l)
    {
      // Write a message about the number of DOFs to be interpolated via
      // surface interpolation.
      if( !(l%writeFreq) )
        std::cout << "Grid zone " << zoneID+1 << ": " << l << " out of "
        << pointsMinDistSearch.size()
        << " DOFs interpolated via minimum distance search."
        << std::endl << std::flush;
      
      // Set a pointer to the coordinates to be searched.
      const su2double *coor = coorInterpol.data() + pointsMinDistSearch[l]*nDim;
      
      // Carry out the minimum distance search,
      unsigned short subElem;
      unsigned long  parElem;
      int            rank;
      su2double      dist;
      su2double      weightsInterpol[4];
      surfaceADT.DetermineNearestElement(coor, dist, subElem, parElem, rank,
                                         weightsInterpol, BBoxTargets,
                                         frontLeaves, frontLeavesNew);
      
      // Subelement found that minimizes the distance to the given coordinate.
      // However, what is needed is the location in the high order parent element.
      // Determine this.
      su2double parCoor[3], wallCoor[3];
      HighOrderMinDistanceSearch(coor, parElem, subElem, weightsInterpol,
                                 &standardBoundaryFacesGrid[indInStandardBoundaryFaces[parElem]],
                                 &surfElems[parElem], coorGrid, parCoor, wallCoor);
      
      // Convert the parametric coordinates of the surface element to
      // parametric weights of the corresponding volume element.
      surfElems[parElem].ConvertParCoorToVolume(&volElems[adjElemID[parElem]],
                                                faceIDInElement[parElem], parCoor);
      
      // Carry out the actual interpolation.
      const unsigned short ind = indInStandardElements[adjElemID[parElem]];
      const unsigned long  ll  = pointsMinDistSearch[l] + zoneOffsetOutputSol;
      
      SolInterpolate(&standardElementsSol[ind], inputSol, zoneOffsetInputSol,
                     &volElems[adjElemID[parElem]], solFormatInput, parCoor,
                     mSolDOFs[ll]);
    }
  
  // Write a message that the surface search is finished.
  std::cout << "Grid zone " << zoneID+1 << ": Surface search finished"
  << std::endl << std::flush;
}

void CFEMInterpolationSol::HighOrderMinDistanceSearch(
                                          const su2double                            *coor,
                                          const unsigned long                        parElem,
                                          const unsigned short                       subElem,
                                          const su2double                            *weightsSubElem,
                                          CFEMStandardBoundaryFace                   *standardBoundaryFaceGrid,
                                          const CFEMInterpolationSurfElem            *surfElem,
                                          const std::vector<std::vector<su2double> > &coorGrid,
                                          su2double                                  *parCoor,
                                          su2double                                  *wallCoor)
{
  // Easier storage of the number of dimensions.
  const unsigned short nDim = coorGrid.size();
  
  // Definition of the maximum number of iterations in the Newton solver
  // and the tolerance level. */
  //const unsigned short maxIt = 50;
  //const su2double tolNewton  = 1.e-10;
  
  /*--------------------------------------------------------------------------*/
  /* Step 1: Create an initial guess for the parametric coordinates from      */
  /*         interpolation in the linear sub-element of the parent element.   */
  /*--------------------------------------------------------------------------*/
  
  // Get the required information for the linear sub-element from the
  // standard boundary face.
  const unsigned short nDOFsPerFace = standardBoundaryFaceGrid->GetNDOFsPerSubFace();
  const unsigned short *connSubFace = standardBoundaryFaceGrid->GetSubFaceConn()
  + subElem*nDOFsPerFace;
  
  // Get the parametric coordinates of the DOFs from the standard boundary face.
  const std::vector<su2double> *locDOFs[] = {standardBoundaryFaceGrid->GetRDOFsFace(),
    standardBoundaryFaceGrid->GetSDOFsFace()};
  
  // Create the initial guess of the parametric coordinates by interpolation
  // in the sub-element. Note that the number of parametric dimensions of a
  // surface element is one less than the number of physical dimensions.
  for(unsigned short iDim=0; iDim<(nDim-1); ++iDim) {
    parCoor[iDim] = 0.0;
    const su2double *coorDOFs = locDOFs[iDim]->data();
    for(unsigned short i=0; i<nDOFsPerFace; ++i)
      parCoor[iDim] += weightsSubElem[i]*coorDOFs[connSubFace[i]];
  }
  
  /*--------------------------------------------------------------------------*/
  /* Step 2: The Newton algorithm to compute the parametric coordinates in    */
  /*         the high order element.                                          */
  /*--------------------------------------------------------------------------*/
  
  // Determine the number of DOFs of the standard element and allocate the
  // memory to store the basis functions and its derivatives.
  const unsigned short nDOFs = standardBoundaryFaceGrid->GetNDOFsFace();
  
  std::vector<su2double> lagBasis(nDOFs);
  std::vector<std::vector<su2double> > dLagBasis(nDim-1, std::vector<su2double>(nDOFs));
  
  // Easier storage of the DOFs of the face.
  const unsigned long *DOFs = surfElem->mConnGrid.data();
  
  // Compute the Lagrangian basis functions and its derivatives in
  // the current parametric coordinate.
  standardBoundaryFaceGrid->FaceBasisFunctionsAndDerivativesInPoint(parCoor, lagBasis,
                                                                    dLagBasis);
  
  // Compute the coordinates of the wall point. Make a distinction between
  // 2D and 3D for efficiency reasons.
  switch( nDim ) {
    case 2: {
      // Two dimensional computation. Compute the values of the coordinates.
      wallCoor[0] = wallCoor[1] = 0.0;
      for(unsigned short i=0; i<nDOFs; ++i) {
        const su2double x = coorGrid[0][DOFs[i]];
        const su2double y = coorGrid[1][DOFs[i]];
        
        wallCoor[0] += x*lagBasis[i];
        wallCoor[1] += y*lagBasis[i];
      }
      
      break;
    }
      
    case 3: {
      // Three dimensional computation. Compute the values of the coordinates.
      wallCoor[0] = wallCoor[1] = wallCoor[2] = 0.0;
      
      for(unsigned short i=0; i<nDOFs; ++i) {
        const su2double x = coorGrid[0][DOFs[i]];
        const su2double y = coorGrid[1][DOFs[i]];
        const su2double z = coorGrid[2][DOFs[i]];
        
        wallCoor[0] += x*lagBasis[i];
        wallCoor[1] += y*lagBasis[i];
        wallCoor[2] += z*lagBasis[i];
      }
      
      break;
    }
  }
  
  // This must still be done.
}

void CFEMInterpolationSol::SolInterpolate(CFEMStandardElement *standardElementSol,
                              const CFEMInterpolationSol      *inputSol,
                              const unsigned long             zoneOffsetInputSol,
                              const CFEMInterpolationVolElem  *volElem,
                              const SolutionFormatT           solFormatInput,
                              const su2double                 *parCoor,
                              std::vector<su2double>          &solDOF)
{
  // Easier storage of the solution DOFs of the input grid.
  const std::vector<std::vector<su2double> > &solInput = inputSol->GetSolDOFs();
  
  // Determine the number of DOFs in the donor element.
  const unsigned short nDOFs = standardElementSol->GetNDOFs();
  
  // Determine the interpolation weights inside this element.
  std::vector<su2double> wSol(nDOFs);
  standardElementSol->BasisFunctionsInPoint(parCoor, wSol);
  
  // Initialize the solution to be interpolated to zero.
  for(unsigned short var=0; var<solDOF.size(); ++var)
    solDOF[var] = 0.0;
  
  // Interpolate the solution, depending on the input solution format.
  switch( solFormatInput )
  {
    case VertexCentered:
    case FEM:
    {
      // Vertex centered or FEM format. Carry out the interpolation.
      for(unsigned short k=0; k<nDOFs; ++k)
      {
        const unsigned long kk = zoneOffsetInputSol + volElem->mConnGrid[k];
        for(unsigned short var=0; var<solDOF.size(); ++var)
          solDOF[var] += wSol[k]*solInput[kk][var];
      }
      
      break;
    }
      
      //--------------------------------------------------------------------------
      
    case CellCentered:
    {
      // Cell centered format. This is a temporary implementation, where
      // the data of the cell center is used.
      const unsigned long kk = zoneOffsetInputSol + volElem->mElemID;
      for(unsigned short var=0; var<solDOF.size(); ++var)
        solDOF[var] = solInput[kk][var];
      
      break;
    }
      
      //--------------------------------------------------------------------------
      
    case DG_FEM:
    {
      // Discontinuous Galerkin format. Carry out the interpolation.
      for(unsigned short k=0; k<nDOFs; ++k)
      {
        const unsigned long kk = zoneOffsetInputSol
        + volElem->mOffsetSolDOFsDG + k;
        for(unsigned short var=0; var<solDOF.size(); ++var)
          solDOF[var] += wSol[k]*solInput[kk][var];
      }
      
      break;
    }
  }
}

CFEMInterpolationFaceOfElem::CFEMInterpolationFaceOfElem()
{
  nCornerPoints   = 0;
  cornerPoints[0] = cornerPoints[1] = cornerPoints[2] = cornerPoints[3] = ULONG_MAX;
  elemID          = ULONG_MAX;
  faceID          = 0;
}

bool CFEMInterpolationFaceOfElem::operator<(const CFEMInterpolationFaceOfElem &other) const
{
  if(nCornerPoints != other.nCornerPoints) return nCornerPoints < other.nCornerPoints;
  
  for(unsigned short i=0; i<nCornerPoints; ++i)
    if(cornerPoints[i] != other.cornerPoints[i]) return cornerPoints[i] < other.cornerPoints[i];
  
  return false;
}

void CFEMInterpolationSurfElem::ConvertParCoorToVolume(const CFEMInterpolationVolElem   *volElem,
                                                       const unsigned short             faceIDInElement,
                                                       su2double                        *parCoor) const
{
  //----------------------------------------------------------------------------
  // Step 1: Transform the parametric coordinates of the face such that it
  //         corresponds to the face numbering used in the volume element.
  //----------------------------------------------------------------------------
  
  // Get the corner points of the surface and the corresponding face in the
  // element.
  unsigned short nCornerPoints;
  unsigned long  cornerPointsSurf[4], cornerPointsVol[4];
  
  GetCornerPoints(nCornerPoints, cornerPointsSurf);
  volElem->GetCornerPointsFace(faceIDInElement, nCornerPoints, cornerPointsVol);
  
  // Determine the element type and set the parametric coordinates r and s
  // accordingly.
  su2double r = 0.0, s = 0.0;
  switch( mVTK_TYPE )
  {
    case LINE:
    {
      // Line element. Either the directions coincide or run opposite.
      if(     (cornerPointsSurf[0] == cornerPointsVol[0]) &&
         (cornerPointsSurf[1] == cornerPointsVol[1])) r = parCoor[0];
      else if((cornerPointsSurf[0] == cornerPointsVol[1]) &&
              (cornerPointsSurf[1] == cornerPointsVol[0])) r = -parCoor[0];
      else
        SU2_MPI::Error("Lines do not match. This should not happen.", CURRENT_FUNCTION);
      break;
    }
      
    case TRIANGLE:
    {
      // Element is a triangle. Determine the orientation and set the parametric
      // coordinates accordingly.
      if(cornerPointsSurf[0] == cornerPointsVol[0])
      {
        if((cornerPointsSurf[1] == cornerPointsVol[1]) &&
           (cornerPointsSurf[2] == cornerPointsVol[2]))
        {
          r = parCoor[0];
          s = parCoor[1];
        }
        else if((cornerPointsSurf[1] == cornerPointsVol[2]) &&
                (cornerPointsSurf[2] == cornerPointsVol[1]))
        {
          r = parCoor[1];
          s = parCoor[0];
        }
        else
          SU2_MPI::Error("Triangles do not match. This should not happen.", CURRENT_FUNCTION);
      }
      else if(cornerPointsSurf[0] == cornerPointsVol[1])
      {
        if((cornerPointsSurf[1] == cornerPointsVol[2]) &&
           (cornerPointsSurf[2] == cornerPointsVol[0]) )
        {
          r = -1.0 - parCoor[0] - parCoor[1];
          s = parCoor[0];
        }
        else if((cornerPointsSurf[1] == cornerPointsVol[0]) &&
                (cornerPointsSurf[2] == cornerPointsVol[2]))
        {
          r = -1.0 - parCoor[0] - parCoor[1];
          s = parCoor[1];
        }
        else
          SU2_MPI::Error("Triangles do not match. This should not happen.", CURRENT_FUNCTION);
      }
      else if(cornerPointsSurf[0] == cornerPointsVol[2])
      {
        if((cornerPointsSurf[1] == cornerPointsVol[0]) &&
           (cornerPointsSurf[2] == cornerPointsVol[1]) )
        {
          r = parCoor[1];
          s = -1.0 - parCoor[0] - parCoor[1];
        }
        else if((cornerPointsSurf[1] == cornerPointsVol[1]) &&
                (cornerPointsSurf[2] == cornerPointsVol[0]))
        {
          r = parCoor[0];
          s = -1.0 - parCoor[0] - parCoor[1];
        }
        else
          SU2_MPI::Error("Triangles do not match. This should not happen.", CURRENT_FUNCTION);
      }
      else
        SU2_MPI::Error("Triangles do not match. This should not happen.", CURRENT_FUNCTION);
      break;
    }
      
    case QUADRILATERAL:
    {
      // Element is a quadrilateral. Determine the orientation and set the parametric
      // coordinates accordingly.
      if(cornerPointsSurf[0] == cornerPointsVol[0])
      {
        if((cornerPointsSurf[1] == cornerPointsVol[1]) &&
           (cornerPointsSurf[2] == cornerPointsVol[2]) &&
           (cornerPointsSurf[3] == cornerPointsVol[3]))
        {
          r = parCoor[0];
          s = parCoor[1];
        }
        else if((cornerPointsSurf[1] == cornerPointsVol[3]) &&
                (cornerPointsSurf[2] == cornerPointsVol[2]) &&
                (cornerPointsSurf[3] == cornerPointsVol[1]))
        {
          r = parCoor[1];
          s = parCoor[0];
        }
        else
          SU2_MPI::Error("Quadrilaterals do not match. This should not happen.", CURRENT_FUNCTION);
      }
      else if(cornerPointsSurf[0] == cornerPointsVol[1])
      {
        if((cornerPointsSurf[1] == cornerPointsVol[0]) &&
           (cornerPointsSurf[2] == cornerPointsVol[3]) &&
           (cornerPointsSurf[3] == cornerPointsVol[2]))
        {
          r = -parCoor[0];
          s =  parCoor[1];
        }
        else if((cornerPointsSurf[1] == cornerPointsVol[2]) &&
                (cornerPointsSurf[2] == cornerPointsVol[3]) &&
                (cornerPointsSurf[3] == cornerPointsVol[0]))
        {
          r = -parCoor[1];
          s =  parCoor[0];
        }
        else
          SU2_MPI::Error("Quadrilaterals do not match. This should not happen.", CURRENT_FUNCTION);
      }
      else if(cornerPointsSurf[0] == cornerPointsVol[2])
      {
        if((cornerPointsSurf[1] == cornerPointsVol[1]) &&
           (cornerPointsSurf[2] == cornerPointsVol[0]) &&
           (cornerPointsSurf[3] == cornerPointsVol[3]))
        {
          r = -parCoor[1];
          s = -parCoor[0];
        }
        else if((cornerPointsSurf[1] == cornerPointsVol[3]) &&
                (cornerPointsSurf[2] == cornerPointsVol[0]) &&
                (cornerPointsSurf[3] == cornerPointsVol[1]))
        {
          r = -parCoor[0];
          s = -parCoor[1];
        }
        else
          SU2_MPI::Error("Quadrilaterals do not match. This should not happen.", CURRENT_FUNCTION);
      }
      else if(cornerPointsSurf[0] == cornerPointsVol[3])
      {
        if((cornerPointsSurf[1] == cornerPointsVol[0]) &&
           (cornerPointsSurf[2] == cornerPointsVol[1]) &&
           (cornerPointsSurf[3] == cornerPointsVol[2]))
        {
          r =  parCoor[1];
          s = -parCoor[0];
        }
        else if((cornerPointsSurf[1] == cornerPointsVol[2]) &&
                (cornerPointsSurf[2] == cornerPointsVol[1]) &&
                (cornerPointsSurf[3] == cornerPointsVol[0]))
        {
          r =  parCoor[0];
          s = -parCoor[1];
        }
      }
      else
        SU2_MPI::Error("Quadrilaterals do not match. This should not happen.", CURRENT_FUNCTION);
      break;
    }
      
    default:
    {
      SU2_MPI::Error("This should not happen.", CURRENT_FUNCTION);
    }
  }
  
  //----------------------------------------------------------------------------
  // Step 2: Transform the parametric face coordinates of the face to
  //         parametric volume coordinates.
  //----------------------------------------------------------------------------
  
  // Determine the volume element type and act accordingly.
  switch( volElem->mVTK_TYPE )
  {
    case TRIANGLE:
    {
      // Volume element is a triangle. Determine the face ID and set the
      // parametric volume coordinates from the surface coordinates.
      switch(faceIDInElement )
      {
        case 0: parCoor[0] =  r;   parCoor[1] = -1.0; break;
        case 1: parCoor[0] = -r;   parCoor[1] =  r;   break;
        case 2: parCoor[0] = -1.0; parCoor[1] = -r;   break;
        default:
          SU2_MPI::Error("Invalid face ID of a triangle. This should not happen.", CURRENT_FUNCTION);
      }
      break;
    }
      
    case QUADRILATERAL:
    {
      // Volume element is a quadrilatral. Determine the face ID and set the
      // parametric volume coordinates from the surface coordinates.
      switch(faceIDInElement )
      {
        case 0: parCoor[0] =  r;   parCoor[1] = -1.0; break;
        case 1: parCoor[0] =  1.0; parCoor[1] =  r;   break;
        case 2: parCoor[0] = -r;   parCoor[1] =  1.0; break;
        case 3: parCoor[0] = -1.0; parCoor[1] = -r;   break;
        default:
          SU2_MPI::Error("Invalid face ID of a quadrilatral. This should not happen.", CURRENT_FUNCTION);
      }
      break;
    }
      
    case TETRAHEDRON:
    {
      // Volume element is a tetrahedron. Determine the face ID and set the
      // parametric volume coordinates from the surface coordinates.
      switch(faceIDInElement )
      {
        case 0: parCoor[0] =  r;       parCoor[1] =  s;   parCoor[2] = -1.0; break;
        case 1: parCoor[0] =  s;       parCoor[1] = -1.0; parCoor[2] =  r;   break;
        case 2: parCoor[0] = -1.0;     parCoor[1] =  r;   parCoor[2] =  s;   break;
        case 3: parCoor[0] = -1.0-r-s; parCoor[1] =  s;   parCoor[2] =  r;   break;
        default:
          SU2_MPI::Error("Invalid face ID of a tetrahedron. This should not happen.", CURRENT_FUNCTION);
      }
      break;
    }
      
    case PYRAMID:
    {
      // Volume element is a pyramid. Determine the face ID and set the
      // parametric volume coordinates from the surface coordinates.
      switch(faceIDInElement )
      {
        case 0: parCoor[0] = r;               parCoor[1] = s;               parCoor[2] = -1.0; break;
        case 1: parCoor[0] = 0.5*(1.0+r) + s; parCoor[1] = 0.5*(r-1.0);     parCoor[2] =  r;   break;
        case 2: parCoor[0] = 0.5*(1.0+s) + r; parCoor[1] = 0.5*(1.0-s);     parCoor[2] =  s;   break;
        case 3: parCoor[0] = 0.5*(s-1.0);     parCoor[1] = 0.5*(1.0+s) + r; parCoor[2] =  s;   break;
        case 4: parCoor[0] = 0.5*(r-1.0);     parCoor[1] = 0.5*(1.0+r) + s; parCoor[2] =  r;   break;
        default:
          SU2_MPI::Error("Invalid face ID of a pyramid. This should not happen.", CURRENT_FUNCTION);
      }
      break;
    }
      
    case PRISM:
    {
      // Volume element is a prism. Determine the face ID and set the
      // parametric volume coordinates from the surface coordinates.
      switch(faceIDInElement )
      {
        case 0: parCoor[0] =  r;   parCoor[1] =  s;   parCoor[2] = -1.0; break;
        case 1: parCoor[0] =  s;   parCoor[1] =  r;   parCoor[2] =  1.0; break;
        case 2: parCoor[0] =  s;   parCoor[1] = -1.0; parCoor[2] =  r;   break;
        case 3: parCoor[0] = -1.0; parCoor[1] =  r;   parCoor[2] =  s;   break;
        case 4: parCoor[0] = -s;   parCoor[1] =  s;   parCoor[2] =  r;   break;
        default:
          SU2_MPI::Error("Invalid face ID of a prism. This should not happen.", CURRENT_FUNCTION);
      }
      break;
    }
      
    case HEXAHEDRON:
    {
      // Volume element is a hexahedron. Determine the face ID and set the
      // parametric volume coordinates from the surface coordinates.
      switch(faceIDInElement )
      {
        case 0: parCoor[0] =  r;   parCoor[1] =  s;   parCoor[2] = -1.0; break;
        case 1: parCoor[0] =  s;   parCoor[1] =  r;   parCoor[2] =  1.0; break;
        case 2: parCoor[0] =  s;   parCoor[1] = -1.0; parCoor[2] =  r;   break;
        case 3: parCoor[0] =  r;   parCoor[1] =  1.0; parCoor[2] =  s;   break;
        case 4: parCoor[0] = -1.0; parCoor[1] =  r;   parCoor[2] =  s;   break;
        case 5: parCoor[0] =  1.0; parCoor[1] =  s;   parCoor[2] =  r;   break;
        default:
          SU2_MPI::Error("Invalid face ID of a hexahedron. This should not happen.", CURRENT_FUNCTION);
      }
      break;
    }
      
    default:
    {
      SU2_MPI::Error("This should not happen.", CURRENT_FUNCTION);
    }
  }
}

void CFEMInterpolationSurfElem::GetCornerPoints(unsigned short &nCornerPoints,
                                                unsigned long  cornerPoints[]) const
{
  // Get the local ID's, relative to the element, of the corner points.
  switch( mVTK_TYPE )
  {
    case LINE:
    {
      nCornerPoints = 2;
      cornerPoints[0] = 0; cornerPoints[1] = mNPolyGrid;
      break;
    }
      
    case TRIANGLE:
    {
      nCornerPoints = 3;
      cornerPoints[0] = 0; cornerPoints[1] = mNPolyGrid; cornerPoints[2] = mNDOFsGrid-1;
      break;
    }
      
    case QUADRILATERAL:
    {
      unsigned short nn2 = mNPolyGrid*(mNPolyGrid+1);
      nCornerPoints = 4;
      cornerPoints[0] = 0;            cornerPoints[1] = mNPolyGrid;
      cornerPoints[2] = mNDOFsGrid-1; cornerPoints[3] = nn2;
      break;
    }
      
    default:
    {
      SU2_MPI::Error("This should not happen.", CURRENT_FUNCTION);
      break;
    }
  }
  
  // Convert the local ID's to global ID's.
  for(unsigned short j=0; j<nCornerPoints; ++j)
    cornerPoints[j] = mConnGrid[cornerPoints[j]];
}

void CFEMInterpolationSurfElem::StoreElemData(const unsigned short VTK_Type,
                                              const unsigned short nPolyGrid,
                                              const unsigned short nDOFsGrid,
                                              const unsigned long  *connGrid)
{
  // Copy the scalar integer data.
  mVTK_TYPE  = VTK_Type;
  mNPolyGrid = nPolyGrid;
  mNDOFsGrid = nDOFsGrid;
  
  // Allocate the memory for the connectivity of the grid.
  mConnGrid.resize(mNDOFsGrid);
  
  // Copy the connectivity data.
  for(unsigned short i=0; i<mNDOFsGrid; ++i)
    mConnGrid[i] = connGrid[i];
  
  // If a linear element is used, the node numbering for non-simplices
  // must be adapted. The reason is that compatability with the original
  // SU2 format is maintained for linear elements, but for the FEM solver
  // the nodes of the elements are stored row-wise.
  if((nPolyGrid == 1) && (VTK_Type == QUADRILATERAL))
    std::swap(mConnGrid[2], mConnGrid[3]);
}

void CFEMInterpolationVolElem::GetCornerPointsAllFaces(unsigned short &nFaces,
                                                       unsigned short nPointsPerFace[],
                                                       unsigned long  faceConn[6][4]) const
{
  // Determine the element type and set the face data accordingly.
  // The faceConn values are local to the element.
  const unsigned short nPoly = mNPolyGrid;
  const unsigned short nDOFs = mNDOFsGrid;
  unsigned short nn2, nn3, nn4;
  switch( mVTK_TYPE )
  {
    case TRIANGLE:
      nFaces = 3;
      nPointsPerFace[0] = 2; faceConn[0][0] = 0;        faceConn[0][1] = nPoly;
      nPointsPerFace[1] = 2; faceConn[1][0] = nPoly;    faceConn[1][1] = nDOFs -1;
      nPointsPerFace[2] = 2; faceConn[2][0] = nDOFs -1; faceConn[2][1] = 0;
      break;
      
    case QUADRILATERAL:
      nFaces = 4; nn2 = nPoly*(nPoly+1);
      nPointsPerFace[0] = 2; faceConn[0][0] = 0;        faceConn[0][1] = nPoly;
      nPointsPerFace[1] = 2; faceConn[1][0] = nPoly;    faceConn[1][1] = nDOFs -1;
      nPointsPerFace[2] = 2; faceConn[2][0] = nDOFs -1; faceConn[2][1] = nn2;
      nPointsPerFace[3] = 2; faceConn[3][0] = nn2;      faceConn[3][1] = 0;
      break;
      
    case TETRAHEDRON:
      nFaces = 4; nn2 = (nPoly+1)*(nPoly+2)/2 -1; nn3 = nDOFs -1;
      nPointsPerFace[0] = 3; faceConn[0][0] = 0;     faceConn[0][1] = nPoly; faceConn[0][2] = nn2;
      nPointsPerFace[1] = 3; faceConn[1][0] = 0;     faceConn[1][1] = nn3;   faceConn[1][2] = nPoly;
      nPointsPerFace[2] = 3; faceConn[2][0] = 0;     faceConn[2][1] = nn2;   faceConn[2][2] = nn3;
      nPointsPerFace[3] = 3; faceConn[3][0] = nPoly; faceConn[3][1] = nn3;   faceConn[3][2] = nn2;
      break;
      
    case PYRAMID:
      nFaces = 5; nn2 = (nPoly+1)*(nPoly+1) -1; nn3 = nn2 - nPoly;
      nPointsPerFace[0] = 4; faceConn[0][0] = 0;     faceConn[0][1] = nPoly;    faceConn[0][2] = nn2; faceConn[0][3] = nn3;
      nPointsPerFace[1] = 3; faceConn[1][0] = 0;     faceConn[1][1] = nDOFs -1; faceConn[1][2] = nPoly;
      nPointsPerFace[2] = 3; faceConn[2][0] = nn3;   faceConn[2][1] = nn2;      faceConn[2][2] = nDOFs -1;
      nPointsPerFace[3] = 3; faceConn[3][0] = 0;     faceConn[3][1] = nn3;      faceConn[3][2] = nDOFs -1;
      nPointsPerFace[4] = 3; faceConn[4][0] = nPoly; faceConn[4][1] = nDOFs -1; faceConn[4][2] = nn2;
      break;
      
    case PRISM:
      nFaces = 5; nn2 = (nPoly+1)*(nPoly+2)/2; nn3 = nPoly*nn2; --nn2;
      nPointsPerFace[0] = 3; faceConn[0][0] = 0;     faceConn[0][1] = nPoly;     faceConn[0][2] = nn2;
      nPointsPerFace[1] = 3; faceConn[1][0] = nn3;   faceConn[1][1] = nn2+nn3;   faceConn[1][2] = nPoly+nn3;
      nPointsPerFace[2] = 4; faceConn[2][0] = 0;     faceConn[2][1] = nn3;       faceConn[2][2] = nPoly+nn3; faceConn[2][3] = nPoly;
      nPointsPerFace[3] = 4; faceConn[3][0] = 0;     faceConn[3][1] = nn2;       faceConn[3][2] = nn2+nn3;   faceConn[3][3] = nn3;
      nPointsPerFace[4] = 4; faceConn[4][0] = nPoly; faceConn[4][1] = nPoly+nn3; faceConn[4][2] = nn2+nn3;   faceConn[4][3] = nn2;
      break;
      
    case HEXAHEDRON:
      nFaces = 6; nn2 = (nPoly+1)*(nPoly+1); nn4 = nPoly*nn2; --nn2; nn3 = nn2 - nPoly;
      nPointsPerFace[0] = 4; faceConn[0][0] = 0;     faceConn[0][1] = nPoly;     faceConn[0][2] = nn2;       faceConn[0][3] = nn3;
      nPointsPerFace[1] = 4; faceConn[1][0] = nn4;   faceConn[1][1] = nn3+nn4;   faceConn[1][2] = nn2+nn4;   faceConn[1][3] = nPoly+nn4;
      nPointsPerFace[2] = 4; faceConn[2][0] = 0;     faceConn[2][1] = nn4;       faceConn[2][2] = nPoly+nn4; faceConn[2][3] = nPoly;
      nPointsPerFace[3] = 4; faceConn[3][0] = nn3;   faceConn[3][1] = nn2;       faceConn[3][2] = nn2+nn4;   faceConn[3][3] = nn3+nn4;
      nPointsPerFace[4] = 4; faceConn[4][0] = 0;     faceConn[4][1] = nn3;       faceConn[4][2] = nn3+nn4;   faceConn[4][3] = nn4;
      nPointsPerFace[5] = 4; faceConn[5][0] = nPoly; faceConn[5][1] = nPoly+nn4; faceConn[5][2] = nn2+nn4;   faceConn[5][3] = nn2;
      break;
      
    default:
      SU2_MPI::Error("This should not happen.", CURRENT_FUNCTION);
      break;
  }
  
  // Convert the local ID's to global ID's.
  for(unsigned short i=0; i<nFaces; ++i)
  {
    for(unsigned short j=0; j<nPointsPerFace[i]; ++j)
      faceConn[i][j] = mConnGrid[faceConn[i][j]];
  }
}

void CFEMInterpolationVolElem::GetCornerPointsFace(const unsigned short faceIDInElement,
                                                   unsigned short       &nCornerPoints,
                                                   unsigned long        cornerPoints[]) const
{
  // This is not the most efficient implementation, but certainly the easiest.
  // Get all the faces of the element.
  unsigned short nFaces;
  unsigned short nPointsPerFace[6];
  unsigned long  faceConn[6][4];
  GetCornerPointsAllFaces(nFaces, nPointsPerFace, faceConn);
  
  // Copy the data for the requested face.
  nCornerPoints = nPointsPerFace[faceIDInElement];
  for(unsigned short i=0; i<nCornerPoints; ++i)
    cornerPoints[i] = faceConn[faceIDInElement][i];
}

void CFEMInterpolationVolElem::StoreElemData(const unsigned long  elemID,
                                             const unsigned short VTK_Type,
                                             const unsigned short nPolyGrid,
                                             const unsigned short nPolySol,
                                             const unsigned short nDOFsGrid,
                                             const unsigned short nDOFsSol,
                                             const unsigned long  offsetSolDOFsDG,
                                             const unsigned long  *connGrid)
{
  // Copy the scalar integer data.
  mElemID          = elemID;
  mVTK_TYPE        = VTK_Type;
  mNPolyGrid       = nPolyGrid;
  mNPolySol        = nPolySol;
  mNDOFsGrid       = nDOFsGrid;
  mNDOFsSol        = nDOFsSol;
  mOffsetSolDOFsDG = offsetSolDOFsDG;
  
  // Allocate the memory for the connectivity of the grid.
  mConnGrid.resize(mNDOFsGrid);
  
  // Copy the connectivity data.
  for(unsigned short i=0; i<mNDOFsGrid; ++i)
    mConnGrid[i] = connGrid[i];
  
  // If a linear element is used, the node numbering for non-simplices
  // must be adapted. The reason is that compatability with the original
  // SU2 format is maintained for linear elements, but for the FEM solver
  // the nodes of the elements are stored row-wise.
  if(nPolyGrid == 1)
  {
    switch( VTK_Type )
    {
      case QUADRILATERAL:
        std::swap(mConnGrid[2], mConnGrid[3]);
        break;
        
      case HEXAHEDRON:
        std::swap(mConnGrid[2], mConnGrid[3]);
        std::swap(mConnGrid[6], mConnGrid[7]);
        break;
        
      case PYRAMID:
        std::swap(mConnGrid[2], mConnGrid[3]);
        break;
    }
  }
}

void CFEMInterpolationGridZone::CopyZoneData(CGeometry* geometry)
{
  unsigned long iElem, nElem, iPoint, nPoint;
  unsigned short iNode, iDim, nDim, iMarker, nMarker;
  
  // Allocate the memory for volume elements.
  nElem = geometry->GetnElem();
  mVolElems.resize(nElem);
  
  // Loop over the elements.
  unsigned long nSolDOFs = 0;
  for(iElem = 0; iElem < nElem; iElem++){
    int VTKType, nPolyGrid, nPolySol, nDOFsGrid, nDOFsSol;
    
    // Determine the VTK type of the element
    VTKType = geometry->elem[iElem]->GetVTK_Type();
    
    // Determine the polynomial degree of the grid and solution. GetNPoly returns 0 if not FEM.
    nPolyGrid = geometry->elem[iElem]->GetNPolyGrid();
    nPolySol  = geometry->elem[iElem]->GetNPolySol();
    if(nPolyGrid == 0){
      nPolyGrid = 1;
      nPolySol  = nPolyGrid;
    }
    
    // Determine the number of DOFs for the grid and solution.
    nDOFsGrid = DetermineNDOFs(VTKType, nPolyGrid);
    nDOFsSol  = DetermineNDOFs(VTKType, nPolySol);
    
    // Allocate the memory for the connectivity and read it.
    std::vector<unsigned long> connSU2(nDOFsGrid);
    for(iNode = 0; iNode < nDOFsGrid; iNode++)
      connSU2[iNode] = geometry->elem[iElem]->GetNode(iNode);
    
    // Store the data for this element.
    mVolElems[iElem].StoreElemData(iElem, VTKType, nPolyGrid, nPolySol, nDOFsGrid,
                                   nDOFsSol, nSolDOFs, connSU2.data());
    
    // Update nSolDOFs.
    nSolDOFs += nDOFsSol;

  }
  
  // Allocate the memory for the coordinates.
  nPoint = geometry->GetnPoint();
  nDim   = geometry->GetnDim();
  mCoor.resize(nDim);
  for(iDim = 0; iDim < nDim; iDim++){
    mCoor[iDim].resize(nPoint);
    // Copy the coordinates.
    for(iPoint = 0; iPoint < nPoint; iPoint++){
      mCoor[iDim][iPoint] = geometry->node[iPoint]->GetCoord()[iDim];
    }
  }
  
  // Allocate the memory for the surface elements.
  nElem = geometry->nLocal_Bound_Elem;
  mSurfElems.resize(nElem);
  
  // Loop over the boundary markers to store the surface connectivity.
  // Note that the surface connectivity is stored as one entity. The
  // information of the boundary markers is not kept.
  nElem = 0;
  nMarker = config->GetnMarker_All();
  for(iMarker = 0; iMarker < nMarker; iMarker++){
    for(iElem = 0 iElem = geometry->GetnElem_Bound(iMarker); iElem++, nElem++){
      int VTKType, nPolyGrid, nDOFsGrid;
      
      // Determine the VTK type of the element
      VTKType = geometry->bound[iMarker][iElem]->GetVTK_Type();
      
      // Determine the polynomial degree of the grid and solution. GetNPoly returns 0 if not FEM.
      nPolyGrid = geometry->bound[iMarker][iElem]->GetNPolyGrid();
      if(nPolyGrid == 0){
        nPolyGrid = 1;
      }
      
      // Determine the number of DOFs for the grid and solution.
      nDOFsGrid = DetermineNDOFs(VTKType, nPolyGrid);
      
      // Allocate the memory for the connectivity and read it.
      std::vector<unsigned long> connSU2(nDOFsGrid);
      for(iNode = 0; iNode < nDOFsGrid; iNode++)
        connSU2[iNode] = geometry->bound[iMarker][iElem]->GetNode(iNode);
      
      // Store the data for this element.
      mSurfElems[nElem].StoreElemData(VTKType, nPolyGrid, nDOFsGrid,
                                      connSU2.data());
    }
  }
  
}

void CFEMInterpolationGridZone::DetermineCoorInterpolation(std::vector<su2double> &coorInterpol,
                                                           const SolutionFormatT  solFormatWrite)
{
  // Determine the number of dimensions.
  const unsigned short nDim = mCoor.size();
  
  // Make a distinction between the requested formats.
  switch( solFormatWrite )
  {
    case VertexCentered:
    case FEM:
    {
      // Vertex centered scheme. The coordinates for interpolation are simply the
      // coordinates of the grid.
      const unsigned long nPoints = mCoor[0].size();
      coorInterpol.resize(nDim*nPoints);
      
      unsigned long ii = 0;
      for(unsigned long i=0; i<nPoints; ++i)
      {
        for(unsigned short iDim=0; iDim<nDim; ++iDim, ++ii)
          coorInterpol[ii] = mCoor[iDim][i];
      }
      
      break;
    }
      
      //--------------------------------------------------------------------------
      
    case CellCentered:
    {
      // Check if high order elements are present.
      if( HighOrderElementsInZone() )
        SU2_MPI::Error("Cell centered data requested, but high order elements present. Not possible.", CURRENT_FUNCTION);
      
      // Allocate the memory for the cell centered coordinates.
      const unsigned long nElem = mVolElems.size();
      coorInterpol.assign(nDim*nElem, 0.0);
      
      // Loop over the elements and compute the coordinates of the cell centers.
      for(unsigned long i=0; i<nElem; ++i)
      {
        su2double *CC = coorInterpol.data() + i*nDim;
        
        for(unsigned short j=0; j<mVolElems[i].mNDOFsGrid; ++j)
        {
          const unsigned long ii = mVolElems[i].mConnGrid[j];
          for(unsigned short iDim=0; iDim<nDim; ++iDim)
            CC[iDim] += mCoor[iDim][ii];
        }
        
        for(unsigned short iDim=0; iDim<nDim; ++iDim)
          CC[iDim] /= mVolElems[i].mNDOFsGrid;
      }
      
      break;
    }
      
      //--------------------------------------------------------------------------
      
    case DG_FEM:
    {
      // Determine the number of DOFs for which the coordinates must be
      // determined.
      const unsigned long nDOFsTot = mVolElems.back().mOffsetSolDOFsDG
      + mVolElems.back().mNDOFsSol;
      
      // Allocate the memory for the coordinates. Initialize it to zero.
      coorInterpol.assign(nDim*nDOFsTot, 0.0);
      
      // Define the vectors of the standard elements in case an
      // interpolation must be carried out for the coordinates.
      std::vector<CFEMStandardElement> standardElementsGrid;
      std::vector<CFEMStandardElement> standardElementsSol;
      
      // Loop over the elements to determine the coordinates of its DOFs.
      unsigned long ii = 0;
      for(unsigned long i=0; i<mVolElems.size(); ++i)
      {
        // Determine if the grid DOFs and solution DOFs coincide.
        if(mVolElems[i].mNDOFsGrid == mVolElems[i].mNDOFsSol)
        {
          // Simply copy the coordinates.
          for(unsigned short j=0; j<mVolElems[i].mNDOFsGrid; ++j)
          {
            const unsigned long jj = mVolElems[i].mConnGrid[j];
            for(unsigned short iDim=0; iDim<nDim; ++iDim, ++ii)
              coorInterpol[ii] = mCoor[iDim][jj];
          }
        }
        else
        {
          // The polynomial degree of the grid and solution differs.
          // The coordinates of the solution DOFs must be interpolated.
          // First determine the index in the list of standard elements.
          // If not present yet, create the standard elements.
          unsigned long j;
          for(j=0; j<standardElementsSol.size(); ++j)
          {
            if(standardElementsSol[j].SameStandardElement(mVolElems[i].mVTK_TYPE,
                                                          mVolElems[i].mNPolySol,
                                                          false) &&
               standardElementsGrid[j].SameStandardElement(mVolElems[i].mVTK_TYPE,
                                                           mVolElems[i].mNPolyGrid,
                                                           false))
              break;
          }
          
          if(j == standardElementsSol.size())
          {
            standardElementsSol.push_back(CFEMStandardElement(mVolElems[i].mVTK_TYPE,
                                                              mVolElems[i].mNPolySol,
                                                              false));
            standardElementsGrid.push_back(CFEMStandardElement(mVolElems[i].mVTK_TYPE,
                                                               mVolElems[i].mNPolyGrid,
                                                               false,
                                                               standardElementsSol[j].GetOrderExact(),
                                                               standardElementsSol[j].GetRDOFs(),
                                                               standardElementsSol[j].GetSDOFs(),
                                                               standardElementsSol[j].GetTDOFs()));
          }
          
          // Get the interpolation data for the coordinates of the solution
          // DOFs. The standard element of the grid must be used to obtain this
          // information.
          const su2double *lagBasisSolDOFs = standardElementsGrid[j].GetBasisFunctionsSolDOFs();
          
          // Loop over the solution DOFs to interpolate the coordinates.
          for(unsigned short j=0; j<mVolElems[i].mNDOFsSol; ++j)
          {
            // Set the pointer to store the interpolation data for this DOF.
            const su2double *weights = lagBasisSolDOFs + j*mVolElems[i].mNDOFsGrid;
            
            // Loop over the number of dimensions.
            for(unsigned short iDim=0; iDim<nDim; ++iDim, ++ii)
            {
              for(unsigned short j=0; j<mVolElems[i].mNDOFsGrid; ++j)
              {
                const unsigned long jj = mVolElems[i].mConnGrid[j];
                coorInterpol[ii] += weights[j]*mCoor[iDim][jj];
              }
            }
          }
        }
      }
      
      break;
    }
  }
}

size_t CFEMInterpolationGridZone::GetNSolDOFs(const SolutionFormatT solFormat) const
{
  // Determine the solution format and return the appropriate value.
  switch( solFormat )
  {
    case FEM:
    case VertexCentered: return GetNGridDOFs();   break;
    case CellCentered:   return GetNVolumeElem(); break;
    case DG_FEM:         return GetNSolDOFsDG();  break;
  }
  
  // Return 0 to avoid a compiler warning.
  return 0;
}

bool CFEMInterpolationGridZone::HighOrderElementsInZone(void) const
{
  // Loop over the volume elements and check if at least one high order
  // element is present.
  bool highOrder = false;
  for(unsigned long i=0; i<mVolElems.size(); ++i)
    if(mVolElems[i].mNPolySol > 1) highOrder = true;
  
  // Return highOrder.
  return highOrder;
}

int CFEMInterpolationGridZone::ReadSU2ZoneData(std::ifstream &su2File,
                                               const int     zoneID,
                                               const bool    multipleZones)
{
  std::ostringstream message;
  
  // Read the spatial dimension for this zone. If no information is found about
  // the number of dimensions, set it to its default value of 3.
  int nDim;
  ResetPositionIfstream(su2File, zoneID, multipleZones);
  if( !FindIntValueFromKeyword(su2File, "ndime", nDim) ) nDim = 3;
  
  // Determine the number of volume elements in this zone.
  int nElem;
  if( !FindIntValueFromKeyword(su2File, "nelem", nElem) )
    SU2_MPI::Error("No volume elements found.", CURRENT_FUNCTION);
  
  // Allocate the memory for volume elements.
  mVolElems.resize(nElem);
  
  // Loop over the elements.
  unsigned long nSolDOFs = 0;
  for(int i=0; i<nElem; ++i)
  {
    // Read the entire line as a string and make it ready for reading.
    std::string lineBuf;
    std::getline(su2File, lineBuf);
    std::istringstream istr(lineBuf);
    
    // Read the su2 element type and determine the polynomial degree and
    // DOFs for both the grid and solution.
    int su2ElemType, VTKType, nPolyGrid, nPolySol, nDOFsGrid, nDOFsSol;
    istr >> su2ElemType;
    
    DetermineElementInfo(su2ElemType, VTKType, nPolyGrid, nPolySol,
                         nDOFsGrid, nDOFsSol);
    
    // Allocate the memory for the connectivity and read it.
    std::vector<unsigned long> connSU2(nDOFsGrid);
    for(int j=0; j<nDOFsGrid; ++j)
      istr >> connSU2[j];
    
    // Store the data for this element.
    mVolElems[i].StoreElemData(i, VTKType, nPolyGrid, nPolySol, nDOFsGrid,
                               nDOFsSol, nSolDOFs, connSU2.data());
    
    // Update nSolDOFs.
    nSolDOFs += nDOFsSol;
  }
  
  // Determine the number of grid points in this zone.
  ResetPositionIfstream(su2File, zoneID, multipleZones);
  int nPoints;
  if( !FindIntValueFromKeyword(su2File, "npoin", nPoints) )
  {
    message << "No points for zone " << zoneID
    << " found in the su2 file.";
    SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }
  
  // Allocate the memory for the coordinates.
  mCoor.resize(nDim);
  for(int iDim=0; iDim<nDim; ++iDim)
    mCoor[iDim].resize(nPoints);
  
  // Read the coordinates.
  for(int i=0; i<nPoints; ++i)
  {
    std::string lineBuf;
    std::getline(su2File, lineBuf);
    std::istringstream istr(lineBuf);
    
    for(int iDim=0; iDim<nDim; ++iDim)
      istr >> mCoor[iDim][i];
  }
  
  // Set the position in ifstream to the beginning of this zone and
  // determine the number of boundary markers.
  ResetPositionIfstream(su2File, zoneID, multipleZones);
  
  int nMarkers;
  if( !FindIntValueFromKeyword(su2File, "nmark", nMarkers) )
  {
    message << "No boundary markers for zone " << zoneID
    << " found in the su2 file.";
    SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }
  
  // Loop over the boundary markers to determine the total number of surface elements.
  nElem = 0;
  for(int iMarker=0; iMarker<nMarkers; ++iMarker)
  {
    // Read the number of surface elements for this marker and update nElem.
    int nSurfElem;
    if( !FindIntValueFromKeyword(su2File, "marker_elems", nSurfElem) )
      SU2_MPI::Error("String \"MARKER_ELEMS\" not found.", CURRENT_FUNCTION);
    nElem += nSurfElem;
  }
  
  // Allocate the memory for the surface elements.
  mSurfElems.resize(nElem);
  
  // Reset the position in ifstream to the location where the marker
  // information for this zone starts.
  ResetPositionIfstream(su2File, zoneID, multipleZones);
  FindIntValueFromKeyword(su2File, "nmark", nMarkers);
  
  // Loop over the boundary markers to store the surface connectivity.
  // Note that the surface connectivity is stored as one entity. The
  // information of the boundary markers is not kept.
  nElem = 0;
  for(int iMarker=0; iMarker<nMarkers; ++iMarker)
  {
    // Read the number of surface elements for this marker.
    int nSurfElem;
    FindIntValueFromKeyword(su2File, "marker_elems", nSurfElem);
    
    // Loop over the surface elements.
    for(int i=0; i<nSurfElem; ++i, ++nElem)
    {
      // Read the entire line as a string and make it ready for reading.
      std::string lineBuf;
      std::getline(su2File, lineBuf);
      std::istringstream istr(lineBuf);
      
      // Read the su2 element type and determine the polynomial degree and
      // DOFs for both the grid and solution. For surface elements the
      // grid and solution information is always the same and therefore
      // nPolySol and nDOFsSol are not used afterwards.
      int su2ElemType, VTKType, nPolyGrid, nPolySol, nDOFsGrid, nDOFsSol;
      istr >> su2ElemType;
      
      DetermineElementInfo(su2ElemType, VTKType, nPolyGrid, nPolySol,
                           nDOFsGrid, nDOFsSol);
      
      // Allocate the memory for the connectivity and read it.
      std::vector<unsigned long> connSU2(nDOFsGrid);
      for(int j=0; j<nDOFsGrid; ++j)
        istr >> connSU2[j];
      
      // Store the data for this surface element.
      mSurfElems[nElem].StoreElemData(VTKType, nPolyGrid, nDOFsGrid,
                                      connSU2.data());
    }
  }
  
  // Return the dimension for this zone.
  return nDim;
}

void CFEMInterpolationGridZone::DetermineElementInfo(int su2ElemType,
                                                     int &VTKType,
                                                     int &nPolyGrid,
                                                     int &nPolySol,
                                                     int &nDOFsGrid,
                                                     int &nDOFsSol)
{
  // Check if the su2ElemType is larger than 10000. If that is the case then
  // the polynomial degree of the grid and solution is different.
  if(su2ElemType > 10000)
  {
    nPolySol    = su2ElemType/10000 - 1;
    su2ElemType = su2ElemType%10000;
    nPolyGrid   = su2ElemType/100 + 1;
  }
  else
  {
    nPolyGrid = su2ElemType/100 + 1;
    nPolySol  = nPolyGrid;
  }
  
  // Determine the VTK type of the element.
  VTKType = su2ElemType%100;
  
  // Determine the number of DOFs for the grid and solution.
  nDOFsGrid = DetermineNDOFs(VTKType, nPolyGrid);
  nDOFsSol  = DetermineNDOFs(VTKType, nPolySol);
}

int CFEMInterpolationGridZone::DetermineNDOFs(const int VTKType,
                                              const int nPoly)
{
  // Initialization.
  const int nDOFsEdge = nPoly + 1;
  int nDOFs = 0;
  
  // Determine the element type and set the number of DOFs from the polynomial
  // degree of the element.
  switch( VTKType )
  {
    case LINE:
      nDOFs = nDOFsEdge;
      break;
      
    case TRIANGLE:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)/2;
      break;
      
    case QUADRILATERAL:
      nDOFs = nDOFsEdge*nDOFsEdge;
      break;
      
    case TETRAHEDRON:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)*(nDOFsEdge+2)/6;
      break;
      
    case PYRAMID:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)*(2*nDOFsEdge+1)/6;
      break;
      
    case PRISM:
      nDOFs = nDOFsEdge*nDOFsEdge*(nDOFsEdge+1)/2;
      break;
      
    case HEXAHEDRON:
      nDOFs = nDOFsEdge*nDOFsEdge*nDOFsEdge;
      break;
      
    default:
      SU2_MPI::Error("Unsupported element type encountered.", CURRENT_FUNCTION);
  }
  
  // Return nDOFs.
  return nDOFs;
}

void CFEMInterpolationGridZone::ResetPositionIfstream(std::ifstream &su2File,
                                                      const int     zoneID,
                                                      const bool    multipleZones)
{
  // Reset the position to the beginning of the file.
  su2File.clear();
  su2File.seekg(0);
  
  // When multiple zones are present, search the correct zone.
  if( multipleZones )
  {
    int zoneNumber = -1;
    while(zoneNumber != zoneID)
    {
      if( !FindIntValueFromKeyword(su2File, "izone", zoneNumber) )
      {
        std::ostringstream message;
        message << "No information for zone " << zoneID
        << " found in the su2 file.";
        SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
      }
    }
  }
}

CFEMInterpolationGrid::CFEMInterpolationGrid(void){}

CFEMInterpolationGrid::CFEMInterpolationGrid(CConfig**      config,
                                             CGeometry**    geometry,
                                             unsigned short nZone)
{
  unsigned short iZone;
  
  // Loop over the number of zones to copy the data from geometry to grid class.
  for(iZone = 0; iZone < nZone; iZone++)
    mGridZones[iZone].CopyZoneData(geometry[iZone]);
}

CFEMInterpolationGrid::~CFEMInterpolationGrid(void){}

void CFEMInterpolationGrid::DetermineCoorInterpolation(std::vector<std::vector<su2double> > &coorInterpol,
                                                       const SolutionFormatT                solFormatWrite)
{
  // Allocate the first index for coorInterpol.
  coorInterpol.resize(mGridZones.size());
  
  // Loop over the zones and determine the coordinates for which the
  // interpolation must be carried out.
  for(unsigned long i=0; i<mGridZones.size(); ++i)
    mGridZones[i].DetermineCoorInterpolation(coorInterpol[i], solFormatWrite);
}

void CFEMInterpolationGrid::DetermineSolutionFormat(const int nSolDOFs)
{
  // Determine the total number of elements, grid DOFs and DG solution DOFs in
  // the grid by looping over the zones in this grid. Also determine
  // whether or not high order elements are present.
  unsigned long nElem = 0, nGridDOFs = 0, nSolDOFsDG = 0;
  bool highOrder = false;
  
  for(std::vector<CFEMInterpolationGridZone>::const_iterator ZI =mGridZones.begin();
      ZI!=mGridZones.end(); ++ZI)
  {
    nElem      += ZI->GetNVolumeElem();
    nGridDOFs  += ZI->GetNGridDOFs();
    nSolDOFsDG += ZI->GetNSolDOFsDG();
    
    if( ZI->HighOrderElementsInZone() ) highOrder = true;
  }
  
  // Determine the format of the solution file.
  if(nSolDOFs == (int) nGridDOFs)
  {
    // A vertex centered scheme is used, i.e. VertexCentered or FEM.
    // When high order elements are present, set the format the FEM,
    // otherwise VertexCentered is most logical (linear FEM gives the
    // same result).
    if( highOrder ) mSolutionFormat = FEM;
    else            mSolutionFormat = VertexCentered;
  }
  else
  {
    // An element based format is used, which is either CellCentered or DG_FEM.
    // For CellCentered the number of solution DOFs is equal to the number of
    // elements. For DG_FEM the number of solution DOFs is equal to nSolDOFsDG.
    if(     nSolDOFs == (int) nElem)      mSolutionFormat = CellCentered;
    else if(nSolDOFs == (int) nSolDOFsDG) mSolutionFormat = DG_FEM;
    else
    {
      // Unknown format. Write an error message and exit.
      SU2_MPI::Error("Unknown solution format.", CURRENT_FUNCTION);
    }
  }
}

