/*
 * Copyright 1997-2003, Regents of the University of Minnesota
 *
 * parmetis.h
 *
 * This file contains function prototypes and constrant definitions for 
 * ParMETIS
 *
 * Started 7/21/03
 * George
 *
 */

#ifndef __parmetis_h__
#define __parmetis_h__

#include <mpi.h>
#include <metis.h>

#ifndef _MSC_VER
#ifndef __MINGW32__
#define __cdecl
#endif
#endif

#if IDXTYPEWIDTH == 32
  /*#define IDX_T         MPI_INT32_T */
  #define IDX_T         MPI_INT
  #define KEEP_BIT      0x40000000L
#elif IDXTYPEWIDTH == 64
  /*#define IDX_T         MPI_INT64_T */
  #define IDX_T         MPI_LONG_LONG_INT
  #define KEEP_BIT      0x4000000000000000LL
#else
  #error "Incorrect user-supplied value fo IDXTYPEWIDTH"
#endif


#if REALTYPEWIDTH == 32
  #define REAL_T        MPI_FLOAT
#elif REALTYPEWIDTH == 64
  #define REAL_T        MPI_DOUBLE
#else
  #error "Incorrect user-supplied value fo REALTYPEWIDTH"
#endif



/*************************************************************************
* Constants 
**************************************************************************/
#define PARMETIS_MAJOR_VERSION        4
#define PARMETIS_MINOR_VERSION        0
#define PARMETIS_SUBMINOR_VERSION     3


/*************************************************************************
* Function prototypes
**************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif

/*-------------------------------------------------------------------
* API Introduced with Release 3.0 (current API) 
*--------------------------------------------------------------------*/
int __cdecl ParMETIS_V3_PartKway(
             idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
	     idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts, 
	     real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part, 
	     MPI_Comm *comm);

int __cdecl ParMETIS_V3_PartGeomKway(
             idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
	     idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ndims, real_t *xyz, 
	     idx_t *ncon, idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, 
	     idx_t *edgecut, idx_t *part, MPI_Comm *comm);

int __cdecl ParMETIS_V3_PartGeom(
             idx_t *vtxdist, idx_t *ndims, real_t *xyz, idx_t *part, MPI_Comm *comm);

int __cdecl ParMETIS_V3_RefineKway(
             idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
	     idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts, 
	     real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, 
	     idx_t *part, MPI_Comm *comm);

int __cdecl ParMETIS_V3_AdaptiveRepart(
             idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
	     idx_t *vsize, idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, 
	     idx_t *nparts, real_t *tpwgts, real_t *ubvec, real_t *ipc2redist, 
	     idx_t *options, idx_t *edgecut, idx_t *part, MPI_Comm *comm);

int __cdecl ParMETIS_V3_Mesh2Dual(
             idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *numflag, 
	     idx_t *ncommonnodes, idx_t **xadj, idx_t **adjncy, MPI_Comm *comm);

int __cdecl ParMETIS_V3_PartMeshKway(
             idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *elmwgt, 
	     idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *ncommonnodes, idx_t *nparts, 
	     real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part, 
	     MPI_Comm *comm);

int __cdecl ParMETIS_V3_NodeND(
             idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *numflag, 
	     idx_t *options, idx_t *order, idx_t *sizes, MPI_Comm *comm);

int __cdecl ParMETIS_V32_NodeND(
             idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
             idx_t *numflag, idx_t *mtype, idx_t *rtype, idx_t *p_nseps, idx_t *s_nseps,
             real_t *ubfrac, idx_t *seed, idx_t *dbglvl, idx_t *order, 
             idx_t *sizes, MPI_Comm *comm);

int __cdecl ParMETIS_SerialNodeND(
             idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *numflag, 
             idx_t *options, idx_t *order, idx_t *sizes, MPI_Comm *comm);

#ifdef __cplusplus
}
#endif


/*------------------------------------------------------------------------
* Enum type definitions 
*-------------------------------------------------------------------------*/
/*! Operation type codes */
typedef enum {
  PARMETIS_OP_KMETIS,
  PARMETIS_OP_GKMETIS,
  PARMETIS_OP_GMETIS,
  PARMETIS_OP_RMETIS,
  PARMETIS_OP_AMETIS,
  PARMETIS_OP_OMETIS,
  PARMETIS_OP_M2DUAL,
  PARMETIS_OP_MKMETIS
} pmoptype_et;


/*************************************************************************
* Various constants used for the different parameters
**************************************************************************/
/* Matching types */
#define PARMETIS_MTYPE_LOCAL     1    /* Restrict matching to within processor vertices */
#define PARMETIS_MTYPE_GLOBAL    2    /* Remote vertices can be matched */

/* Separator refinement types */
#define PARMETIS_SRTYPE_GREEDY    1    /* Vertices are visted from highest to lowest gain */
#define PARMETIS_SRTYPE_2PHASE    2    /* Separators are refined in a two-phase fashion using
                                          PARMETIS_SRTYPE_GREEDY for the 2nd phase */

/* Coupling types for ParMETIS_V3_RefineKway & ParMETIS_V3_AdaptiveRepart */
#define PARMETIS_PSR_COUPLED    1    /* # of partitions == # of processors */
#define PARMETIS_PSR_UNCOUPLED  2    /* # of partitions != # of processors */


/* Debug levels (fields should be ORed) */
#define PARMETIS_DBGLVL_TIME        1      /* Perform timing analysis */
#define PARMETIS_DBGLVL_INFO        2      /* Perform timing analysis */
#define PARMETIS_DBGLVL_PROGRESS    4      /* Show the coarsening progress */
#define PARMETIS_DBGLVL_REFINEINFO  8      /* Show info on communication during folding */
#define PARMETIS_DBGLVL_MATCHINFO   16     /* Show info on matching */
#define PARMETIS_DBGLVL_RMOVEINFO   32     /* Show info on communication during folding */
#define PARMETIS_DBGLVL_REMAP       64     /* Determines if remapping will take place */

#endif 
