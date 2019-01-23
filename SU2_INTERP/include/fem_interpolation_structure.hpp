/*!
 * \file fem_interpolation_structure.hpp
 * \brief Headers of the main subroutines for creating the interpolation structure for the FEM solver.
 *        The subroutines and functions are in the <i>fem_interpolation_structure.cpp</i> file.
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

#pragma once

#include "../../Common/include/adt_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/fem_geometry_structure.hpp"
#include "../../SU2_CFD/include/solver_structure.hpp"

using namespace std;

/*!
 * \brief Enumerated type for the solution format.
 */
enum SolutionFormatT
{
  VertexCentered = 1,
  CellCentered   = 2,
  FEM            = 3,
  DG_FEM         = 4
};

/*!
 * \class CFEMInterpolationVolElem
 * \brief Class to store volume elements for interpolation of FEM solution.
 * \version 6.1.0 "Falcon"
 */class CFEMInterpolationVolElem {
public:
  unsigned long mElemID;                 /*!< \brief ID of this element. */
  unsigned short mVTK_TYPE;              /*!< \brief Element type using the VTK convention. */
  unsigned short mNPolyGrid;             /*!< \brief Polynomial degree for the geometry of the element. */
  unsigned short mNPolySol;              /*!< \brief Polynomial degree for the solution of the element. */
  unsigned short mNDOFsGrid;             /*!< \brief Number of DOFs for the geometry of the element. */
  unsigned short mNDOFsSol;              /*!< \brief Number of DOFs for the solution of the element. */
  unsigned long mOffsetSolDOFsDG;        /*!< \brief Offset for the solution DOFs, only for DG formulation. */
  std::vector<unsigned long> mConnGrid;  /*!< \brief The node numbers for the grid DOFs. */
  /*!
   * \brief Constructor of the class.
   */
   CFEMInterpolationVolElem(void){}
  
  /*!
   * \brief Destructor of the class.
   */
   ~CFEMInterpolationVolElem(void){}
  
  /*!
   * \brief Copy constructor.
   */
  CFEMInterpolationVolElem(const CFEMInterpolationVolElem &other) {Copy(other);}
  
  /*!
   * \brief Assignment operator.
   */
  CFEMInterpolationVolElem& operator=(const CFEMInterpolationVolElem &other) {Copy(other); return (*this);}
  
  /*!
   * \brief Function, which makes available the corner points of all faces of the element.
   */
  void GetCornerPointsAllFaces(unsigned short &nFaces,
                               unsigned short nPointsPerFace[],
                               unsigned long  faceConn[6][4]) const;
  
  /*!
   * \brief Function, which makes available the corner points of the requested face of the element.
   */
  void GetCornerPointsFace(const unsigned short faceIDInElement,
                           unsigned short       &nCornerPoints,
                           unsigned long        cornerPoints[]) const;
  
  /*!
   * \brief Function, which sets the data of this element.
   */
  void StoreElemData(const unsigned long  elemID,
                     const unsigned short VTK_Type,
                     const unsigned short nPolyGrid,
                     const unsigned short nPolySol,
                     const unsigned short nDOFsGrid,
                     const unsigned short nDOFsSol,
                     const unsigned long  offsetSolDOFsDG,
                     const unsigned long  *connGrid);
  
private:
  /*!
   * \brief Function, which makes a deep copy.
   */
  void Copy(const CFEMInterpolationVolElem &other);
};

/*!
 * \class CFEMInterpolationFaceOfElem
 * \brief Class to store face of element for interpolation of FEM solution.
 * \version 6.1.0 "Falcon"
 */
class CFEMInterpolationFaceOfElem {
public:
  unsigned short nCornerPoints;          /*!< \brief Number of corner points of the face. */
  unsigned long  cornerPoints[4];        /*!< \brief ID's of the corner points. */
  unsigned long  elemID;                 /*!< \brief Element ID. */
  unsigned short faceID;                 /*!< \brief Face ID inside the element ID. */
  
  /*!
   * \brief Constructor of the class.
   */
  CFEMInterpolationFaceOfElem(void);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CFEMInterpolationFaceOfElem(void){}
  
  /*!
   * \brief Copy constructor.
   */
  CFEMInterpolationFaceOfElem(const CFEMInterpolationFaceOfElem &other){Copy(other);}
  
  /*!
   * \brief Assignment operator.
   */
  CFEMInterpolationFaceOfElem& operator=(const CFEMInterpolationFaceOfElem &other){Copy(other); return (*this);}
  
  /*!
   * \brief Less than operator. Needed for the sorting and searching.
   */
  bool operator<(const CFEMInterpolationFaceOfElem &other) const;
  
  /*!
   * \brief Member function, which creates a unique numbering for the corner points. A sort in increasing order is OK for this purpose.
   */
  void CreateUniqueNumbering(void){std::sort(cornerPoints, cornerPoints+nCornerPoints);}
  
private:
  /*!
   * \brief Copy function, which copies the data of the given object into the current object.
   */
  void Copy(const CFEMInterpolationFaceOfElem &other);
};

/*!
 * \class CFEMInterpolationSurfElem
 * \brief Class to store surface elements for interpolation of FEM solution.
 * \version 6.1.0 "Falcon"
 */
class CFEMInterpolationSurfElem {
public:
  unsigned short mVTK_TYPE;              /*!< \brief Element type using the VTK convention. */
  unsigned short mNPolyGrid;             /*!< \brief Polynomial degree for the geometry of the element. */
  unsigned short mNDOFsGrid;             /*!< \brief Number of DOFs for the geometry of the element. */
  std::vector<unsigned long> mConnGrid;  /*!< \brief The node numbers for the grid DOFs. */
  
  /*!
   * \brief Constructor of the class.
   */
  CFEMInterpolationSurfElem(void){}
  
  /*!
   * \brief Destructor of the class.
   */
  ~CFEMInterpolationSurfElem(void){}
  
  /*!
   * \brief Copy constructor.
   */
  CFEMInterpolationSurfElem(const CFEMInterpolationSurfElem &other) {Copy(other);}
  
  /*!
   * \brief Assignment operator.
   */
  CFEMInterpolationSurfElem& operator=(const CFEMInterpolationSurfElem &other) {Copy(other); return (*this);}
  
  /*!
   * \brief Function, which converts the parametric coordinates of the surface element to the parametric weights of the adjacent volume element.
   */
  void ConvertParCoorToVolume(const CFEMInterpolationVolElem   *volElem,
                              const unsigned short             faceIDInElement,
                              su2double                        *parCoor) const;
  
  /*!
   * \brief Function, which makes available the corner points of the surface element.
   */
  void GetCornerPoints(unsigned short &nCornerPoints,
                       unsigned long  cornerPoints[]) const;
  
  /*!
   * \brief Function, which sets the data of this element.
   */
  void StoreElemData(const unsigned short VTK_Type,
                     const unsigned short nPolyGrid,
                     const unsigned short nDOFsGrid,
                     const unsigned long  *connGrid);
  
private:
  /*!
   * \brief Function, which makes a deep copy.
   */
  void Copy(const CFEMInterpolationSurfElem &other);
};

/*!
 * \class CFEMInterpolationGridZone
 * \brief Class to store grid zone information for interpolation of FEM solution.
 * \version 6.1.0 "Falcon"
 */
class CFEMInterpolationGridZone {
public:
  std::vector<CFEMInterpolationSurfElem> mSurfElems;  /*!< \brief Vector of the surface elements. */
  std::vector<CFEMInterpolationVolElem> mVolElems;    /*!< \brief Vector of the volume elements. */
  std::vector<std::vector<su2double> > mCoor;         /*!< \brief Vector of vector of the coordinates of the DOFs.
                                                       * The first index is the number of spatial dimensions. */
  /*!
   * \brief Constructor of the class.
   */
  CFEMInterpolationGridZone(void){}
  
  /*!
   * \brief Destructor of the class.
   */
  ~CFEMInterpolationGridZone(void){}
  
  /*!
   * \brief Copy constructor.
   */
  CFEMInterpolationGridZone(const CFEMInterpolationGridZone &other) {Copy(other);}
  
  /*!
   * \brief Assignment operator.
   */
  CFEMInterpolationGridZone& operator=(const CFEMInterpolationGridZone &other) {Copy(other); return (*this);}
  
  /*!
   * \brief Function, which copies data from the SU2 geometry structure to the interpolation zone data structure for a single zone.
   */
  void CopySU2GeometryToGrid(CConfig*   config,
                             CGeometry* geometry);
  
  /*!
   * \brief Function, which determines the coordinates for the points to be interpolated.
   */
  void DetermineCoorInterpolation(CConfig*               config,
                                  std::vector<su2double> &coorInterpol,
                                  const SolutionFormatT  solFormatWrite);
  
  /*!
   * \brief Function, which makes available the number of grid DOFs.
   */
  size_t GetNGridDOFs(void) const {return mCoor[0].size();}
  
  /*!
   * \brief Function, which determines the number of solution DOFs, depending on the given solution format.
   */
  size_t GetNSolDOFs(const SolutionFormatT solFormat) const;
  
  /*!
   * \brief Function, which makes available the number of DG solution DOFs.
   */
  size_t GetNSolDOFsDG(void) const {return mVolElems.back().mOffsetSolDOFsDG + mVolElems.back().mNDOFsSol;}
  
  /*!
   * \brief Function, which makes available the number of volume elements.
   */
  size_t GetNVolumeElem(void) const {return mVolElems.size();}
  
  /*!
   * \brief Function, which determines whether or not high order elements are present in the zone.
   */
  bool HighOrderElementsInZone(void) const;
  
private:
  /*!
   * \brief Function, which makes a deep copy.
   */
  void Copy(const CFEMInterpolationGridZone &other);
  
  /*!
   * \brief Function, which determines the information for the given SU2 element type.
   */
  void DetermineElementInfo(int su2ElemType,
                            int &VTKType,
                            int &nPolyGrid,
                            int &nPolySol,
                            int &nDOFsGrid,
                            int &nDOFsSol);
  
  /*!
   * \brief Function, which determines the number of DOFs for the given element type and polynomial degree.
   */
  int DetermineNDOFs(const int VTKType,
                     const int nPoly);
  
};

/*!
 * \class CFEMInterpolationGrid
 * \brief Class to store grid for interpolation of the FEM solution.
 * \version 6.1.0 "Falcon"
 */
class CFEMInterpolationGrid {
public:
  /*!
   * \brief Constructor of the class.
   */
  CFEMInterpolationGrid(void);
  
  /*!
   * \overload
   * \brief Constructor of the class.
   * \param[in] config   - Definition of the particular problem.
   * \param[in] geometry - The grid for interpolation.
   */
  CFEMInterpolationGrid(CConfig**      config,
                        CGeometry***    geometry,
                        unsigned short nZone);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CFEMInterpolationGrid(void);
  
  /*!
   * \brief Function, which determines the coordinates for the points to be interpolated.
   */
  void DetermineCoorInterpolation(CConfig**                            config,
                                  std::vector<std::vector<su2double> > &coorInterpol,
                                  const SolutionFormatT                solFormatWrite);
  
  /*!
   * \brief Function, which determines the solution format from the grid information and the given number of solution DOFs.
   */
  void DetermineSolutionFormat(const int nSolDOFs);
  
  /*!
   * \brief Function, which makes available the number of dimensions.
   */
  unsigned short GetNDim(void) const {return mNDim;}
  
  /*!
   * \brief Function, which makes available the number of zones.
   */
  size_t GetNZones(void) const {return mGridZones.size();}
  
  /*!
   * \brief Make available the given zone as a const pointer.
   */
  const CFEMInterpolationGridZone *GetGridZone(const unsigned short zone) const {return &mGridZones[zone];}
  
  /*!
   * \brief Make available the solution format of this grid.
   */
  SolutionFormatT GetSolutionFormat(void) const {return mSolutionFormat;}
  
private:
  unsigned short mNDim;                               /*!< \brief Number of spatial dimensions. */
  std::vector<CFEMInterpolationGridZone> mGridZones;  /*!< \brief Vector of grid zones. */
  SolutionFormatT mSolutionFormat;                    /*!< \brief The solution format of the corresponding solution. */
};

/*!
 * \class CFEMInterpolationSol
 * \brief Class to perform interpolation of the FEM solution.
 * \version 6.1.0 "Falcon"
 */
class CFEMInterpolationSol {
public:
  /*!
   * \brief Constructor of the class.
   */
  CFEMInterpolationSol(void);
  
  /*!
   * \overload
   * \brief Constructor of the class.
   * \param[in] config          - Definition of the particular problem.
   * \param[in] input_geometry  - The source grid for interpolation.
   * \param[in] output_geometry - The target grid for interpolation.
   * \param[in] input_solution  - The source solution structure for interpolation.
   * \param[in] output_solution - The target solution structure for interpolation.
   */
  CFEMInterpolationSol(CConfig**      config,
                       CGeometry***   input_geometry,
                       CGeometry***    output_geometry,
                       CSolver****      input_solution,
                       CSolver****     output_solution,
                       unsigned short nZone);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CFEMInterpolationSol(void);
  
  /*!
   * \brief Make available the number of variables in the solution file.
   * \return  Number of variables in the solution file.
   */
  size_t GetNVar(void) const {return mVarNames.size();}
  
  /*!
   * \brief Main function for the interpolation of the given coordinates in the given grid and solution.
   */
  void InterpolateSolution(CConfig**                                  config,
                           const std::vector<std::vector<su2double> > &coorInterpol,
                           const CFEMInterpolationGrid                *inputGrid,
                           const CFEMInterpolationSol                 *inputSol,
                           const CFEMInterpolationGrid                *outputGrid);
  
  /*!
   * \brief Make available the solution DOFs.
   * \return  Number of DOFs in the solution.
   */
  const std::vector<std::vector<su2double> > &GetSolDOFs(void) const {return mSolDOFs;}
    
private:
  int mHeaderFile[5];                             /*!< \brief The 5 integers of the header of the solution file. */
  int mIterNumber;                                /*!< \brief The iteration number stored in the solution file. */
  double mMetaData[8];                            /*!< \brief The meta data stored in the solution file. */
  std::vector<std::string> mVarNames;             /*!< \brief The names of the solution variables. */
  std::vector<std::vector<su2double> > mSolDOFs;  /*!< \brief Double vector, which contains the solution variables in the DOFs. */
  
  /*!
   * \brief Function, which applies a curvature correction to the coordinates to obtain a better representation of these coordinates on the target grid.
   */
  void ApplyCurvatureCorrection(CConfig*                         config,
                                const unsigned short             zoneID,
                                const unsigned short             nDim,
                                const CFEMInterpolationGridZone  *inputGridZone,
                                const CFEMInterpolationGridZone  *outputGridZone,
                                const std::vector<su2double>     &coorOriginal,
                                std::vector<su2double>           &coorCorrected);
  
  /*!
   * \brief Function, which builds the ADT of a surface grid.
   */
  void BuildSurfaceADT(CConfig*                                  config,
                       const CFEMInterpolationGridZone           *gridZone,
                       CADTElemClass                             &surfaceADT,
                       std::vector<CFEMStandardBoundaryFace>     &standardBoundaryFacesGrid,
                       std::vector<CFEMStandardBoundaryFace>     &standardBoundaryFacesSol,
                       std::vector<unsigned short>               &indInStandardBoundaryFaces,
                       std::vector<unsigned long>                &adjElemID,
                       std::vector<unsigned short>               &faceIDInElement);
  
  /*!
   * \brief Function, which performs the containment search in a high order element.
   */
  void HighOrderContainmentSearch(const su2double                            *coor,
                                  const unsigned long                        parElem,
                                  const unsigned short                       subElem,
                                  const su2double                            *weightsSubElem,
                                  CFEMStandardElement                        *standardElementGrid,
                                  const CFEMInterpolationVolElem             *volElem,
                                  const std::vector<std::vector<su2double> > &coorGrid,
                                  su2double                                  *parCoor);
  
  /*!
   * \brief Function, which performs the minimum distance search for a high order surface element.
   */
  void HighOrderMinDistanceSearch(const su2double                            *coor,
                                  const unsigned long                        parElem,
                                  const unsigned short                       subElem,
                                  const su2double                            *weightsSubElem,
                                  CFEMStandardBoundaryFace                   *standardBoundaryFaceGrid,
                                  const CFEMInterpolationSurfElem            *surfElem,
                                  const std::vector<std::vector<su2double> > &coorGrid,
                                  su2double                                  *parCoor,
                                  su2double                                  *wallCoor);
  
  /*!
   * \brief Function, which carries out the actual interpolation of the solution.
   */
  void SolInterpolate(CFEMStandardElement             *standardElementSol,
                      const CFEMInterpolationSol      *inputSol,
                      const unsigned long             zoneOffsetInputSol,
                      const CFEMInterpolationVolElem  *volElem,
                      const SolutionFormatT           solFormatInput,
                      const su2double                 *parCoor,
                      std::vector<su2double>          &solDOF);
  
  /*!
   * \brief Function, which carries out the surface minimum distance search for the interpolation.
   */
  void SurfaceInterpolationSolution(CConfig*                             config,
                                    const unsigned short                 zoneID,
                                    const std::vector<su2double>         &coorInterpol,
                                    const CFEMInterpolationGridZone      *gridZone,
                                    const CFEMInterpolationSol           *inputSol,
                                    const unsigned long                  zoneOffsetInputSol,
                                    const unsigned long                  zoneOffsetOutputSol,
                                    const SolutionFormatT                solFormatInput,
                                    const std::vector<unsigned long>     &pointsMinDistSearch,
                                    std::vector<CFEMStandardElement>     &standardElementsSol,
                                    const std::vector<unsigned short>    &indInStandardElements);
  
  /*!
   * \brief Function, which carries out the volume containment search for the interpolation.
   */
  void VolumeInterpolationSolution(CConfig*                             config,
                                   const unsigned short                 zoneID,
                                   const std::vector<su2double>         &coorInterpol,
                                   const CFEMInterpolationGridZone      *gridZone,
                                   const CFEMInterpolationSol           *inputSol,
                                   const unsigned long                  zoneOffsetInputSol,
                                   const unsigned long                  zoneOffsetOutputSol,
                                   const SolutionFormatT                solFormatInput,
                                   std::vector<unsigned long>           &pointsSearchFailed,
                                   std::vector<CFEMStandardElement>     &standardElementsGrid,
                                   std::vector<CFEMStandardElement>     &standardElementsSol,
                                   std::vector<unsigned short>          &indInStandardElements);
  
  /*!
   * \brief Function, which copies data from the interpolation sol data structure to the SU2 solution structure.
   */
  void CopySolToSU2Solution(CConfig**      config,
                            CGeometry***    geometry,
                            CSolver****     solution,
                            unsigned short nZone);
};
