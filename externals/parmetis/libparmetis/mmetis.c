/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mmetis.c
 *
 * This is the entry point of ParMETIS_V3_PartMeshKway
 *
 * Started 10/19/96
 * George
 *
 * $Id: mmetis.c 10573 2011-07-14 13:31:54Z karypis $
 *
 */

#include <parmetislib.h>


/***********************************************************************************
* This function is the entry point of the parallel k-way multilevel mesh partitionioner. 
* This function assumes nothing about the mesh distribution.
* It is the general case.
************************************************************************************/
int ParMETIS_V3_PartMeshKway(idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *elmwgt, 
        idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *ncommon, idx_t *nparts, 
	real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part, 
	MPI_Comm *comm)
{
  idx_t i, status, nvtxs, nedges, gnedges, npes, mype;
  idx_t *xadj, *adjncy;
  ctrl_t *ctrl;
  size_t curmem;

  /* Check the input parameters and return if an error */
  status = CheckInputsPartMeshKway(elmdist, eptr, eind, elmwgt, wgtflag, numflag,
               ncon, ncommon, nparts, tpwgts, ubvec, options, edgecut, part, comm);
  if (GlobalSEMinComm(*comm, status) == 0)
    return METIS_ERROR;

  status = METIS_OK;
  gk_malloc_init();
  curmem = gk_GetCurMemoryUsed();

  /* Setup the ctrl */
  ctrl = SetupCtrl(PARMETIS_OP_MKMETIS, NULL, 1, 1, NULL, NULL, *comm);
  npes = ctrl->npes;
  mype = ctrl->mype;


  /* Create the dual graph */
  STARTTIMER(ctrl, ctrl->MoveTmr);

  ParMETIS_V3_Mesh2Dual(elmdist, eptr, eind, numflag, ncommon, &xadj, &adjncy, 
      &(ctrl->comm));

  if (ctrl->dbglvl&DBG_INFO) {
    nvtxs = elmdist[mype+1]-elmdist[mype];
    nedges = xadj[nvtxs] + (*numflag == 0 ? 0 : -1);
    rprintf(ctrl, "Completed Dual Graph -- Nvtxs: %"PRIDX", Nedges: %"PRIDX" \n", 
            elmdist[npes], GlobalSESum(ctrl, nedges));
  }

  STOPTIMER(ctrl, ctrl->MoveTmr);


  /* Partition the dual graph */
  STARTTIMER(ctrl, ctrl->TotalTmr);

  status = ParMETIS_V3_PartKway(elmdist, xadj, adjncy, elmwgt, NULL, wgtflag, 
               numflag, ncon, nparts, tpwgts, ubvec, options, edgecut, part, 
               &(ctrl->comm));

  STOPTIMER(ctrl, ctrl->TotalTmr);

  IFSET(ctrl->dbglvl, DBG_TIME, PrintTimer(ctrl, ctrl->MoveTmr,	 "   Mesh2Dual"));
  IFSET(ctrl->dbglvl, DBG_TIME, PrintTimer(ctrl, ctrl->TotalTmr, "    ParMETIS"));

  METIS_Free(xadj);
  METIS_Free(adjncy);

  FreeCtrl(&ctrl);
  if (gk_GetCurMemoryUsed() - curmem > 0) {
    printf("ParMETIS appears to have a memory leak of %zdbytes. Report this.\n",
        (ssize_t)(gk_GetCurMemoryUsed() - curmem));
  }
  gk_malloc_cleanup(0);

  return (int)status;
}
