/*!
 * \file gelite.cpp
 * \brief Functions used to call the GELite library. This library is made
 *        available by Pointwise to carry out a projection of a point onto
 *        the geometry definition.
 * \author E. van der Weide
 * \version 4.1.2 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#ifdef HAVE_GELITE

#include "../include/gelite.hpp"


Int32 BuildBSPTree (ProjectionDatabase *const projectionDatabase)
/*==============================================================================
  ABSTRACT
    This routine builds a BSP tree from a ProjectionBSPTree.
  
  PARAMETERS
    projectionDatabase  INOUT
      ProjectionBSPTree
      
    returns:
      0 - no errors occurred
      2 - unspecified failure
  
  REVISIONS
    02/09/05 - J. P. - Modified to use ProjectionBSPTree classes
    04/18/05 - J. P. - Modified to use EntityList for all lists of entities
    05/16/07 - J. P. - Removed check for NULL returned from operator new
    02/08/13 - J. P. - Modified as follows:
                       - Modified ProjectionBSPTree::Build_BSPTree() to take
                         0 for each of the input values, meaning that an appropriate
                         value will be computed
                       - Bug #747
    02/25/13 - J. P. - TopologyProjectionBSPTreeWrapper:
                       - Consolidated interfaces for adding Geometry, Curved
                         Topology, and Mesh Topology to a ProjectionBSPTree
                       - Bug #1241
    
  AUTHOR
    J. P. Abelanet - 08/16/04
==============================================================================*/
{
  Real64        secs;
  Time          time;
  Int32         i;

  time = Time ();

  /* Create the ProjectionBSPTree */
  projectionDatabase->projectionBSPTree = new ProjectionBSPTree;

  /* Add the entities to the database */
  for (i = 0; i < projectionDatabase->entities.Size (); i++)
    if (CurvedModel::Downcast (projectionDatabase->entities[i])) {
      if (projectionDatabase->FaceCoordinates) {
        EntityList<CurvedFace>  Faces;
        Int32           j;

        CurvedModel::Downcast (projectionDatabase->entities[i])->Inquire_Faces (&Faces);
        for (j = 0; j < Faces.Size (); j++)
          if (TopologyProjectionBSPTreeWrapper::Add_Entity (projectionDatabase->projectionBSPTree,
                    Faces[j]) != Error::No_errors)
            return (2);
      } else if (TopologyProjectionBSPTreeWrapper::Add_Entity (projectionDatabase->projectionBSPTree,
            projectionDatabase->entities[i]) != Error::No_errors)
        return (2);
    } else {
      if (TopologyProjectionBSPTreeWrapper::Add_Entity (projectionDatabase->projectionBSPTree,
            projectionDatabase->entities[i]) != Error::No_errors)
        return (2);
    }

  secs = time.ElapsedSecs ();
  ErrorLogger::LogInfo (CharString::Sprintf (
                "Load Time (secs)              = %.2f\n", secs));
  time = Time ();

  /* Build the BSP tree */
  if (projectionDatabase->projectionBSPTree->Build_BSPTree (
                0 /* MaxLevel */, 0 /*MaxElPerNode*/) != Error::No_errors)
    return (2);

  // Output tree statistics
  secs = time.ElapsedSecs ();
  ErrorLogger::LogInfo (CharString::Sprintf (
                "Build Time (secs)             = %.2f\n", secs));
  projectionDatabase->projectionBSPTree->Print_TreeStatistics ();

  return (0);
}

Real64 Entity_Inquire_ProjectionTolerance (const IsectProjPointEnd *entity)
/*==============================================================================
  ABSTRACT
    This routine finds the expected accuracy of a projection onto the given
    entity.
  
  PARAMETERS
    entity  IN
      entity to use
      
    returns:
      Expected accuracy of projections
  
  REVISIONS
    02/02/06 - J. P. - Splitting Topology classes into Curved and Mesh subclasses

  AUTHOR
    J. P. Abelanet - 08/15/04
==============================================================================*/
#define DEFAULT_POINT_PROJECTION_ACCURACY       (Tolerance::GetSamePoint ()/10.0)       /* same as in SurfaceProjectPoint::Compute_CoordMinimumDistance */
{
  Real64        Tolerance;

  if (entity->subEntity != NULL) {
    if (CurvedFace::Downcast (entity->subEntity))
      return (DEFAULT_POINT_PROJECTION_ACCURACY);
    else if (CurvedEdge::Downcast (entity->subEntity)) {
      CurvedEdge::Downcast (entity->subEntity)->Inquire_Tolerance (&Tolerance);
      return (Tolerance);
    } else if (CurvedVertex::Downcast (entity->subEntity)) {
      CurvedVertex::Downcast (entity->subEntity)->Inquire_Tolerance (&Tolerance);
      return (Tolerance);
    } else {
      ErrorLogger::LogError (CharString::Sprintf (
                "ERROR:  Entity_Inquire_ProjectionTolerance:  Unexpected subentity (%s)!\n",
                entity->subEntity->Inquire_Description ()));
      return (0.0);
    }
  } else if (Surface::Downcast (entity->entity))
    return (DEFAULT_POINT_PROJECTION_ACCURACY);
  else if (Curve::Downcast (entity->entity))
    return (DEFAULT_POINT_PROJECTION_ACCURACY);
  else if (CurvedFace::Downcast (entity->entity))
    return (DEFAULT_POINT_PROJECTION_ACCURACY);
  else if (CurvedEdge::Downcast (entity->entity)) {
    CurvedEdge::Downcast (entity->entity)->Inquire_Tolerance (&Tolerance);
    return (Tolerance);
  } else if (Point::Downcast (entity->entity))
    return (0.0);
  else {
    ErrorLogger::LogError (CharString::Sprintf (
                "ERROR:  Entity_Inquire_ProjectionTolerance:  Unexpected entity (%s)!\n",
        entity->entity->Inquire_Description ()));
    return (0.0);
  }
}

void OutputProjectionStats(const IsectProjPoint *projection,
                           const IsectProjPointEnd *BasePoint,
                           const bool error)
/*==============================================================================
  ABSTRACT
    This routine outputs statistics about a projection.
  
  PARAMETERS
    projection  IN
      point that was projected onto from a point relative to BasePoint
      
    BasePoint  IN
      Base point that was offset -
      if NULL, no initial point was offset
      
    error  IN
      true if an error occurred
  
  REVISIONS
    07/12/11 - J. P. - Modified to check for valid Entity -
                       Bug #983
    11/18/14 - J. P. - Modified to use new ToCharString() methods -
                       Bug #1820
                       
  AUTHOR
    J. P. Abelanet - 10/06/03
==============================================================================*/
{
  CharString    charString;

  /* If an error occurred, report that now */
  if (error)
    charString += "ERROR:\n";

  /* Report the base point that was offset */
  if (BasePoint != NULL) {
    charString += "Base:\n";
    charString += BasePoint->ToCharString();
  }

  /* Report the point that was projected onto from the offset point */
  charString += "Proj:\n";
  if (projection->End1.entity != NULL)
    charString += projection->End1.ToCharString();
  else
    charString += "none\n";

  /* Report the offset point that was projected, as well as the distances */
  if (BasePoint != NULL) {
    charString += "Offset:\n";
    charString += CharString::Sprintf (
                "\tDistBO = %f\n\tDistOP = %f\n\tDistBP = %f\n",
                BasePoint->P.Distance (projection->End2.P),
                projection->Distance,
                projection->End1.entity != NULL?
                        BasePoint->P.Distance (projection->End1.P):
                        Tolerance::Infinity);
  } else {
    charString += "Coord:\n";
    charString += CharString::Sprintf (
                "\tDistance = %f\n",
                projection->Distance);
  }

  charString += "\tP = ";
  charString += projection->End2.P.ToCharString ();
  charString += "\n\n";

  if (error)
    ErrorLogger::LogError (charString);
  else
    ErrorLogger::LogInfo (charString);
}

Int32 ProjectOneCoordinate(Vector3D& Coordinate,
                           ProjectionDatabase *const projectionDatabase)
/*==============================================================================
  ABSTRACT
    This routine projects one point onto the ProjectionBSPTree.
  
  PARAMETERS
    Coordinate  INOUT
      On input:  Coordinate to project
      On output: Projected coordinate.
      
    projectionDatabase  INOUT
      ProjectionBSPTree
  
    returns:
      0 - no errors occurred
      1 - incorrect projection computed
      2 - unspecified projection failure
  
  REVISIONS
    02/09/05 - J. P. - Modified to use ProjectionBSPTree classes
    07/12/11 - J. P. - Modified parameter list to OutputProjectionStats() -
                       Bug #983
    02/25/14 - J. P. - Cleaning up casts involving Entities -
                       Bug #7
    
  AUTHOR
    J. P. Abelanet - 08/16/04
==============================================================================*/
{
  IsectProjPoint        projection;
  Real64                ProjectionTol;
  Error                 error;
  Int32                 ErrorCode = 0;

  /* Project Coordinate onto the database */
  if (projectionDatabase->projectionBSPTree)
    error = projectionDatabase->projectionBSPTree->Compute_CoordMinimumDistance (Coordinate, NULL, NULL,
                &projection);
  else if (projectionDatabase->SingleEntity->Is_InClassID (Geometry::Static_ClassID ()))
    error = Geometry::Downcast (projectionDatabase->SingleEntity)->Compute_CoordMinimumDistance (Coordinate,
                NULL, NULL, &projection);
  else {
    ErrorLogger::LogError (CharString::Sprintf (
                "ERROR:  ProjectOneCoordinate:  Unexpected entity (%s)!\n",
        projectionDatabase->SingleEntity->Inquire_Description ()));
    return (2);
  }

  /* Get the expected accuracy of the projection */
  if (projection.End1.entity != NULL)
    ProjectionTol = Entity_Inquire_ProjectionTolerance (&projection.End1);
  else
    ProjectionTol = 0.0;

  /* Did an error occur? */
  if (error != Error::No_errors)
    ErrorCode = 2;
  else if (projection.Distance > FABS (projectionDatabase->MaxOffset) + ProjectionTol)
    ErrorCode = 1;

  /* Report what happened */
  //if (ErrorCode != 0 || projectionDatabase->Verbose)
  //  OutputProjectionStats (&projection, NULL, ErrorCode != 0);

  // Copy the projected coordinate into Coordinate.
  Coordinate[Vector3D::iX] = projection.End1.P[Vector3D::iX];
  Coordinate[Vector3D::iY] = projection.End1.P[Vector3D::iY];
  Coordinate[Vector3D::iZ] = projection.End1.P[Vector3D::iZ];

  return (ErrorCode);
}

Int32 ReadGeomDefFile (const char *FileName, ProjectionDatabase *const projectionDatabase)
/*==============================================================================
  ABSTRACT
    This reads a file and populates the ProjectionBSPTree.
  
  PARAMETERS
    FileName  IN
      Name of file to read
      
    projectionDatabase  INOUT
      ProjectionBSPTree
      
    returns:
      0 - no errors occurred
      -1 - could not find usable entities
  
  REVISIONS
    05/16/07 - J. P. - Removed check for NULL returned from operator new
    12/12/08 - J. P. - Modified to ignore Planes and to call
                       ForeignEntityPostprocessor::Database_Detach_EntitiesWithDependents() -
                       Bug #307

  AUTHOR
    J. P. Abelanet - 08/16/04
==============================================================================*/
{
  Int32 i;

  /* Create an empty Database */
  projectionDatabase->database = new Database;

  /* Read the file */
  if (NativeTopologyReader::Read (FileName, projectionDatabase->database) != Error::No_errors)
    return (-1);

  // Remove all Entities that have dependent children from the Database,
  //   leaving only the top-level Entities
  ForeignEntityPostprocessor::Database_Detach_EntitiesWithDependents (projectionDatabase->database);

  /* Load all Geometry and Topology other than infinite Planes */
  EntityList<Entity>    entities;
  projectionDatabase->database->Inquire_Entities (entities);
  for (i = 0; i < entities.Size (); i++)
    if (entities[i]->Is_InClassID (Geometry::Static_ClassID ()) ||
        entities[i]->Is_InClassID (Topology::Static_ClassID ()))
      if (!entities[i]->Is_InClassID (Plane::Static_ClassID ()))
        projectionDatabase->entities.Append (entities[i]);

  /* Did we find any usable entities? */
  if (projectionDatabase->entities.Size () == 0)
    return (-1);

  return (0);
}

#endif
