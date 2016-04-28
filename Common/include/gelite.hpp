/*!
 * \file gelite.hpp
 * \brief Include files and headers of the functions used to call the GELite
 *        library. This library is made available by Pointwise to carry out
 *        a projection of a point onto the geometry definition. The functions
 *        are in the <i>gelite.cpp</i> file.
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

#pragma once

#ifdef HAVE_GELITE

/*----------------------------------------------------------------------*/
/*---                    GELite include files.                       ---*/
/*----------------------------------------------------------------------*/

#include <geom/BSplineCurve.h>
#include <geom/BSplineSurface.h>
#include <geom/Curve.h>
#include <geom/Database.h>
#include <geom/Entity.h>
#include <geom/EntityList.h>
#include <geom/ErrorLogger.h>
#include <geom/Geometry.h>
#include <geom/GeometrySession.h>
#include <geom/IsectProjPoint.h>
#include <geom/NativeFileWriter.h>
#include <geom/Plane.h>
#include <geom/Point.h>
#include <geom/ProjectionBSPTree.h>
#include <geom/Surface.h>
#include <nmb/CurvedEdge.h>
#include <nmb/CurvedFace.h>
#include <nmb/CurvedModel.h>
#include <nmb/CurvedVertex.h>
#include <nmb/ForeignEntityPostprocessor.h>
#include <nmb/MeshModel.h>
#include <nmb/NativeTopologyReader.h>
#include <nmb/Topology.h>
#include <nmb/TopologyProjectionBSPTreeWrapper.h>
#include <vmath/Base.h>
#include <vmath/CharString.h>
#include <vmath/Error.h>
#include <vmath/Interval.h>
#include <vmath/IntervalVector2D.h>
#include <vmath/Time.h>
#include <vmath/Tolerance.h>
#include <vmath/Vector2D.h>
#include <vmath/Vector3D.h>

/*!
 * \class ProjectionDatabase
 * \brief Class used for the projection of a point onto the data base. It is
 *        copied from the examples provided with the GELite library.
 * \author J. P. Abelanet (developer of GELite)
 * \version 4.1.2 "Cardinal"
 */
using namespace GE;
class ProjectionDatabase
{
public:
  //------------------------------------------
  // Public constructor.
  //------------------------------------------

  ProjectionDatabase ():
        database (NULL),
        entities (),
        SingleEntity (NULL),
        projectionBSPTree (NULL),
        NumProjections (0),
        NumOutOfTol (0),
        TotalNumProjections (0),
        TotalNumOutOfTol (0),
        MaxOffset (0.0),
        SurfaceSamples (0),
        CurveSamples (0),
        FaceCoordinates (false),
        Verbose (false),
        CreatePoints (false),
        diagnosticEntities ()   {}

  //------------------------------------------
  // Public data members.
  //------------------------------------------

  Database              *database;              /* Database of entities                                 */
  EntityList<Entity>    entities;               /* list of top-level entities being projected onto      */

  Entity                *SingleEntity;          /* single entity to project onto                */
  ProjectionBSPTree     *projectionBSPTree;     /* ProjectionBSPTree - only used if non-NULL    */

  Int32                 NumProjections;         /* current entity - number of projections                       */
  Int32                 NumOutOfTol;            /* current entity - number of projections out of tolerance      */

  Int32                 TotalNumProjections;    /* total number of projections                  */
  Int32                 TotalNumOutOfTol;       /* total number of projections out of tolerance */

  //------------------------------------------
  // User specified options.
  //------------------------------------------

  Real64                MaxOffset;              /* Maximum offset of points relative to the evaluted point -
                                                   all points will be offset a distance in range
                                                   [-MaxOffset, MaxOffset]                      */
  Int32                 SurfaceSamples;         /* # of samples across each surface in U/V      */
  Int32                 CurveSamples;           /* # of samples across each curve in S          */
  bool          FaceCoordinates;        /* if true, convert to face coordinates for topology projections        */
  bool          Verbose;                /* if true, report on each projection           */
  bool          CreatePoints;           /* if true, create a point for each projection  */
  EntityList<Entity>    diagnosticEntities;     /* Diagnostic Entities                          */
};

/*-----------------------------------------------------------------------------*/
/*--- Prototypes of the functions used to interact with the GELite library. ---*/
/*-----------------------------------------------------------------------------*/

Int32 BuildBSPTree(ProjectionDatabase *const projectionDatabase);

Real64 Entity_Inquire_ProjectionTolerance(const IsectProjPointEnd *entity);

void OutputProjectionStats(const IsectProjPoint *projection,
                           const IsectProjPointEnd *BasePoint,
                           const bool error);

Int32 ProjectOneCoordinate(Vector3D& Coordinate,
                           ProjectionDatabase *const projectionDatabase);

Int32 ReadGeomDefFile(const char         *FileName,
                      ProjectionDatabase *const projectionDatabase);

#endif
