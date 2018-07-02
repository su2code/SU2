/*
 * frename.c
 *
 * This file contains some renaming routines to deal with different
 * Fortran compilers.
 *
 * Started 6/1/98
 * George
 *
 * $Id: frename.c 13945 2013-03-30 14:38:24Z karypis $
 *
 */

#include <parmetislib.h>


/*************************************************************************
* Renaming macro (at least to save some typing :))  
**************************************************************************/
#define FRENAME(name0, name1, name2, name3, name4, dargs, cargs)   \
  int name1 dargs {MPI_Comm comm = MPI_Comm_f2c(*icomm); return name0 cargs; }\
  int name2 dargs {MPI_Comm comm = MPI_Comm_f2c(*icomm); return name0 cargs; }\
  int name3 dargs {MPI_Comm comm = MPI_Comm_f2c(*icomm); return name0 cargs; }\
  int name4 dargs {MPI_Comm comm = MPI_Comm_f2c(*icomm); return name0 cargs; }




/*************************************************************************
* Renames for Release 3.0 API
**************************************************************************/
FRENAME(ParMETIS_V3_AdaptiveRepart, 
        PARMETIS_V3_ADAPTIVEREPART,
	parmetis_v3_adaptiverepart,
	parmetis_v3_adaptiverepart_,
	parmetis_v3_adaptiverepart__,
	(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
	 idx_t *vsize, idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon,
	 idx_t *nparts, real_t *tpwgts, real_t *ubvec, real_t *ipc2redist,
	 idx_t *options, idx_t *edgecut, idx_t *part, MPI_Fint *icomm),
	(vtxdist, xadj, adjncy, vwgt, vsize, adjwgt, wgtflag, numflag, ncon,
	 nparts, tpwgts, ubvec, ipc2redist, options, edgecut, part, &comm)
)

FRENAME(ParMETIS_V3_PartGeomKway,
        PARMETIS_V3_PARTGEOMKWAY,
	parmetis_v3_partgeomkway,
	parmetis_v3_partgeomkway_,
	parmetis_v3_partgeomkway__,
        (idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
	 idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ndims, real_t *xyz, 
	 idx_t *ncon, idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, 
	 idx_t *edgecut, idx_t *part, MPI_Fint *icomm),
        (vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, ndims, xyz, 
	 ncon, nparts, tpwgts, ubvec, options, edgecut, part, &comm)
)
	 
FRENAME(ParMETIS_V3_PartGeom,
        PARMETIS_V3_PARTGEOM,
	parmetis_v3_partgeom,
	parmetis_v3_partgeom_,
	parmetis_v3_partgeom__,
	(idx_t *vtxdist, idx_t *ndims, real_t *xyz, idx_t *part, MPI_Fint *icomm),
	(vtxdist, ndims, xyz, part, &comm)
)

FRENAME(ParMETIS_V3_PartKway,
        PARMETIS_V3_PARTKWAY,
	parmetis_v3_partkway,
	parmetis_v3_partkway_,
	parmetis_v3_partkway__,
	(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t *adjwgt, 
	 idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts, real_t *tpwgts, real_t *ubvec, 
	 idx_t *options, idx_t *edgecut, idx_t *part, MPI_Fint *icomm),
	(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, ncon, nparts, tpwgts, 
	 ubvec, options, edgecut, part, &comm)
)

FRENAME(ParMETIS_V3_Mesh2Dual,
        PARMETIS_V3_MESH2DUAL,
	parmetis_v3_mesh2dual,
	parmetis_v3_mesh2dual_,
	parmetis_v3_mesh2dual__,
	(idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *numflag, idx_t *ncommonnodes, 
	 idx_t **xadj, idx_t **adjncy, MPI_Fint *icomm),
	(elmdist, eptr, eind, numflag, ncommonnodes, xadj, adjncy, &comm)
)

FRENAME(ParMETIS_V3_PartMeshKway,
        PARMETIS_V3_PARTMESHKWAY, 
	parmetis_v3_partmeshkway,
	parmetis_v3_partmeshkway_,
	parmetis_v3_partmeshkway__,
	(idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *elmwgt, idx_t *wgtflag, 
	 idx_t *numflag, idx_t *ncon, idx_t *ncommonnodes, idx_t *nparts, real_t *tpwgts, 
	 real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part, MPI_Fint *icomm),
	(elmdist, eptr, eind, elmwgt, wgtflag, numflag, ncon, ncommonnodes, nparts, tpwgts, 
	 ubvec, options, edgecut, part, &comm)
)
	 
FRENAME(ParMETIS_V3_NodeND,
        PARMETIS_V3_NODEND,
        parmetis_v3_nodend,
        parmetis_v3_nodend_,
        parmetis_v3_nodend__,
        (idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *numflag, idx_t *options, 
	 idx_t *order, idx_t *sizes, MPI_Fint *icomm),
        (vtxdist, xadj, adjncy, numflag, options, order, sizes, &comm)
)

FRENAME(ParMETIS_V3_RefineKway,
        PARMETIS_V3_REFINEKWAY,
        parmetis_v3_refinekway,
        parmetis_v3_refinekway_,
        parmetis_v3_refinekway__,
        (idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t *adjwgt, 
	 idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts, real_t *tpwgts, real_t *ubvec, 
	 idx_t *options, idx_t *edgecut, idx_t *part, MPI_Fint *icomm),
        (vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, ncon, nparts, tpwgts, 
	 ubvec, options, edgecut, part, &comm)
)

