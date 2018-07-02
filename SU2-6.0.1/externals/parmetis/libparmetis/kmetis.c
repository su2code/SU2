/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * kmetis.c
 *
 * This is the entry point of ParMETIS_PartKway
 *
 * Started 10/19/96
 * George
 *
 * $Id: kmetis.c 10757 2011-09-15 22:07:47Z karypis $
 *
 */

#include <parmetislib.h>

/***********************************************************************************
* This function is the entry point of the parallel k-way multilevel partitionioner. 
* This function assumes nothing about the graph distribution.
* It is the general case.
************************************************************************************/
int ParMETIS_V3_PartKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
        idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts, 
	real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part, 
	MPI_Comm *comm)
{
  idx_t h, i, status, nvtxs, npes, mype, seed, dbglvl;
  ctrl_t *ctrl;
  graph_t *graph;
  idx_t moptions[METIS_NOPTIONS];
  size_t curmem;


  /* Check the input parameters and return if an error */
  status = CheckInputsPartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, 
               numflag, ncon, nparts, tpwgts, ubvec, options, edgecut, part, 
               comm);
  if (GlobalSEMinComm(*comm, status) == 0) 
    return METIS_ERROR;

  status = METIS_OK;
  gk_malloc_init();
  curmem = gk_GetCurMemoryUsed();

  /* Set up the control */
  ctrl = SetupCtrl(PARMETIS_OP_KMETIS, options, *ncon, *nparts, tpwgts, ubvec, *comm);
  npes = ctrl->npes;
  mype = ctrl->mype;


  /* Take care the nparts == 1 case */
  if (*nparts == 1) {
    iset(vtxdist[mype+1]-vtxdist[mype], (*numflag == 0 ? 0 : 1), part); 
    *edgecut = 0;
    goto DONE;
  }


  /* Take care of npes == 1 case */
  if (npes == 1) {
    nvtxs = vtxdist[1] - vtxdist[0];
    METIS_SetDefaultOptions(moptions);
    moptions[METIS_OPTION_NUMBERING] = *numflag;

    status = METIS_PartGraphKway(&nvtxs, ncon, xadj, adjncy, vwgt, NULL, 
                 adjwgt, nparts, tpwgts, ubvec, moptions, edgecut, part);
 
    goto DONE;
  }


  /* Setup the graph */
  if (*numflag > 0) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  graph = SetupGraph(ctrl, *ncon, vtxdist, xadj, vwgt, NULL, adjncy, adjwgt, *wgtflag);

  /* Setup the workspace */
  AllocateWSpace(ctrl, 10*graph->nvtxs);


  /* Partition the graph */
  STARTTIMER(ctrl, ctrl->TotalTmr);

  ctrl->CoarsenTo = gk_min(vtxdist[npes]+1, 25*(*ncon)*gk_max(npes, *nparts));
  if (vtxdist[npes] < SMALLGRAPH 
      || vtxdist[npes] < npes*20 
      || GlobalSESum(ctrl, graph->nedges) == 0) { /* serially */
    IFSET(ctrl->dbglvl, DBG_INFO, 
        rprintf(ctrl, "Partitioning a graph of size %"PRIDX" serially\n", vtxdist[npes]));
    PartitionSmallGraph(ctrl, graph);
  }
  else { /* in parallel */
    Global_Partition(ctrl, graph);
  }
  ParallelReMapGraph(ctrl, graph);

  icopy(graph->nvtxs, graph->where, part);
  *edgecut = graph->mincut;

  STOPTIMER(ctrl, ctrl->TotalTmr);


  /* Print out stats */
  IFSET(ctrl->dbglvl, DBG_TIME, PrintTimingInfo(ctrl));
  IFSET(ctrl->dbglvl, DBG_TIME, gkMPI_Barrier(ctrl->gcomm));
  IFSET(ctrl->dbglvl, DBG_INFO, PrintPostPartInfo(ctrl, graph, 0));

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



/*************************************************************************/
/*! This function is the driver to the multi-constraint partitioning 
    algorithm.
*/
/*************************************************************************/
void Global_Partition(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, ncon, nparts;
  real_t ftmp, ubavg, lbavg, *lbvec;

  WCOREPUSH;
 
  ncon   = graph->ncon;
  nparts = ctrl->nparts;
  ubavg  = ravg(graph->ncon, ctrl->ubvec);

  CommSetup(ctrl, graph);

  lbvec = rwspacemalloc(ctrl, ncon);

  if (ctrl->dbglvl&DBG_PROGRESS) {
    rprintf(ctrl, "[%6"PRIDX" %8"PRIDX" %5"PRIDX" %5"PRIDX"] [%"PRIDX"] [", graph->gnvtxs, GlobalSESum(ctrl, graph->nedges),
	    GlobalSEMin(ctrl, graph->nvtxs), GlobalSEMax(ctrl, graph->nvtxs), ctrl->CoarsenTo);
    for (i=0; i<ncon; i++)
      rprintf(ctrl, " %.3"PRREAL"", GlobalSEMinFloat(ctrl,graph->nvwgt[rargmin_strd(graph->nvtxs, graph->nvwgt+i, ncon)*ncon+i]));  
    rprintf(ctrl, "] [");
    for (i=0; i<ncon; i++)
      rprintf(ctrl, " %.3"PRREAL"", GlobalSEMaxFloat(ctrl, graph->nvwgt[rargmax_strd(graph->nvtxs, graph->nvwgt+i, ncon)*ncon+i]));  
    rprintf(ctrl, "]\n");
  }

  if (graph->gnvtxs < 1.3*ctrl->CoarsenTo ||
	(graph->finer != NULL &&
	graph->gnvtxs > graph->finer->gnvtxs*COARSEN_FRACTION)) {

    /* Done with coarsening. Find a partition */
    AllocateRefinementWorkSpace(ctrl, 2*graph->nedges);
    graph->where = imalloc(graph->nvtxs+graph->nrecv, "graph->where");

    InitPartition(ctrl, graph);

    if (ctrl->dbglvl&DBG_PROGRESS) {
      ComputePartitionParams(ctrl, graph);
      ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
      rprintf(ctrl, "nvtxs: %10"PRIDX", cut: %8"PRIDX", balance: ", 
          graph->gnvtxs, graph->mincut);
      for (i=0; i<graph->ncon; i++) 
        rprintf(ctrl, "%.3"PRREAL" ", lbvec[i]);
      rprintf(ctrl, "\n");

      /* free memory allocated by ComputePartitionParams */
      gk_free((void **)&graph->ckrinfo, &graph->lnpwgts, &graph->gnpwgts, LTERM);
    }

    /* In case no coarsening took place */
    if (graph->finer == NULL) {
      ComputePartitionParams(ctrl, graph);
      KWayFM(ctrl, graph, NGR_PASSES);
    }
  }
  else {
    Match_Global(ctrl, graph);

    Global_Partition(ctrl, graph->coarser);

    ProjectPartition(ctrl, graph);

    ComputePartitionParams(ctrl, graph);

    if (graph->ncon > 1 && graph->level < 3) {
      for (i=0; i<ncon; i++) {
        ftmp = rsum(nparts, graph->gnpwgts+i, ncon);
        if (ftmp != 0.0)
          lbvec[i] = (real_t)(nparts) *
          graph->gnpwgts[rargmax_strd(nparts, graph->gnpwgts+i, ncon)*ncon+i]/ftmp;
        else
          lbvec[i] = 1.0;
      }
      lbavg = ravg(graph->ncon, lbvec);

      if (lbavg > ubavg + 0.035) {
        if (ctrl->dbglvl&DBG_PROGRESS) {
          ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
          rprintf(ctrl, "nvtxs: %10"PRIDX", cut: %8"PRIDX", balance: ", 
              graph->gnvtxs, graph->mincut);
          for (i=0; i<graph->ncon; i++) 
            rprintf(ctrl, "%.3"PRREAL" ", lbvec[i]);
          rprintf(ctrl, " [b]\n");
	}

        KWayBalance(ctrl, graph, graph->ncon);
      }
    }

    KWayFM(ctrl, graph, NGR_PASSES);

    if (ctrl->dbglvl&DBG_PROGRESS) {
      ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
      rprintf(ctrl, "nvtxs: %10"PRIDX", cut: %8"PRIDX", balance: ", 
          graph->gnvtxs, graph->mincut);
      for (i=0; i<graph->ncon; i++) 
        rprintf(ctrl, "%.3"PRREAL" ", lbvec[i]);
      rprintf(ctrl, "\n");
    }

    if (graph->level != 0)
      gk_free((void **)&graph->lnpwgts, (void **)&graph->gnpwgts, LTERM);
  }

  WCOREPOP;

}


