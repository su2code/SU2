/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * rmetis.c
 *
 * This is the entry point of the partitioning refinement routine
 *
 * Started 10/19/96
 * George
 *
 * $Id: rmetis.c 10612 2011-07-21 20:09:05Z karypis $
 *
 */

#include <parmetislib.h>



/***********************************************************************************
* This function is the entry point of the parallel multilevel local diffusion
* algorithm. It uses parallel undirected diffusion followed by adaptive k-way 
* refinement. This function utilizes local coarsening.
************************************************************************************/
int ParMETIS_V3_RefineKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
        idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts, 
        real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part, 
        MPI_Comm *comm)
{
  idx_t npes, mype, status;
  ctrl_t *ctrl=NULL;
  graph_t *graph=NULL;
  size_t curmem;


  /* Check the input parameters and return if an error */
  status = CheckInputsPartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, 
               numflag, ncon, nparts, tpwgts, ubvec, options, edgecut, part, comm);
  if (GlobalSEMinComm(*comm, status) == 0) 
    return METIS_ERROR;

  status = METIS_OK;
  gk_malloc_init();
  curmem = gk_GetCurMemoryUsed();

  /* Setup ctrl */
  ctrl = SetupCtrl(PARMETIS_OP_RMETIS, options, *ncon, *nparts, tpwgts, ubvec, *comm);
  npes = ctrl->npes;
  mype = ctrl->mype;


  /* Take care the nparts == 1 case */
  if (*nparts == 1) {
    iset(vtxdist[mype+1]-vtxdist[mype], (*numflag == 0 ? 0 : 1), part); 
    *edgecut = 0;
    goto DONE;
  }


  /* setup the graph */
  if (*numflag > 0) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  graph = SetupGraph(ctrl, *ncon, vtxdist, xadj, vwgt, NULL, adjncy, adjwgt, *wgtflag);

  if (ctrl->ps_relation == PARMETIS_PSR_COUPLED)
    iset(graph->nvtxs, mype, graph->home);
  else
    icopy(graph->nvtxs, part, graph->home);


  /* Allocate workspace */
  AllocateWSpace(ctrl, 10*graph->nvtxs);


  /* Partition and Remap */
  STARTTIMER(ctrl, ctrl->TotalTmr);

  ctrl->CoarsenTo = gk_min(vtxdist[npes]+1, 50*(*ncon)*gk_max(npes, *nparts));

  Adaptive_Partition(ctrl, graph);
  ParallelReMapGraph(ctrl, graph);

  icopy(graph->nvtxs, graph->where, part);
  *edgecut = graph->mincut;

  STOPTIMER(ctrl, ctrl->TotalTmr);

  /* Take care of output */
  IFSET(ctrl->dbglvl, DBG_TIME, PrintTimingInfo(ctrl));
  IFSET(ctrl->dbglvl, DBG_TIME, gkMPI_Barrier(ctrl->gcomm));
  IFSET(ctrl->dbglvl, DBG_INFO, PrintPostPartInfo(ctrl, graph, 1));

  FreeInitialGraphAndRemap(graph);

  if (*numflag > 0)
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

DONE:
  FreeCtrl(&ctrl);
  if (gk_GetCurMemoryUsed() - curmem > 0) {
    printf("ParMETIS appears to have a memory leak of %zdbytes. Report this.\n",
        (ssize_t)(gk_GetCurMemoryUsed() - curmem));
  }
  gk_malloc_cleanup(0);

  return (int)status;
}


