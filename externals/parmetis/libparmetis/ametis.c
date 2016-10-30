/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * ametis.c
 *
 * This is the entry point of parallel difussive repartitioning routines
 *
 * Started 10/19/96
 * George
 *
 * $Id: ametis.c 10757 2011-09-15 22:07:47Z karypis $
 *
 */

#include <parmetislib.h>


/***********************************************************************************
* This function is the entry point of the parallel multilevel local diffusion
* algorithm. It uses parallel undirected diffusion followed by adaptive k-way 
* refinement. This function utilizes local coarsening.
************************************************************************************/
int ParMETIS_V3_AdaptiveRepart(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy,
        idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag,
        idx_t *ncon, idx_t *nparts, real_t *tpwgts, real_t *ubvec, real_t *ipc2redist,
        idx_t *options, idx_t *edgecut, idx_t *part, MPI_Comm *comm)
{
  idx_t i, npes, mype, status;
  ctrl_t *ctrl=NULL;
  graph_t *graph=NULL;
  size_t curmem;


  /* Check the input parameters and return if an error */
  status = CheckInputsAdaptiveRepart(vtxdist, xadj, adjncy, vwgt, vsize, adjwgt,
               wgtflag, numflag, ncon, nparts, tpwgts, ubvec, ipc2redist, options,
               edgecut, part, comm);
  if (GlobalSEMinComm(*comm, status) == 0) 
    return METIS_ERROR;

  status = METIS_OK;
  gk_malloc_init();
  curmem = gk_GetCurMemoryUsed();

  /* Setup the ctrl */
  ctrl = SetupCtrl(PARMETIS_OP_AMETIS, options, *ncon, *nparts, tpwgts, ubvec, *comm);
  npes = ctrl->npes;
  mype = ctrl->mype;


  /* Take care the nparts == 1 case */
  if (*nparts == 1) {
    iset(vtxdist[mype+1]-vtxdist[mype], (*numflag == 0 ? 0 : 1), part); 
    *edgecut = 0;
    goto DONE;
  }


  /* Setup the graph */
  if (*numflag > 0) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  graph = SetupGraph(ctrl, *ncon, vtxdist, xadj, vwgt, vsize, adjncy, adjwgt, *wgtflag);

  if (ctrl->ps_relation == PARMETIS_PSR_COUPLED)
    iset(graph->nvtxs, mype, graph->home);
  else {
    /* Downgrade the partition numbers if part[] has more partitions that nparts */
    for (i=0; i<graph->nvtxs; i++)
      part[i] = (part[i] >= ctrl->nparts ? 0 : part[i]);

    icopy(graph->nvtxs, part, graph->home);
  }


  /* Allocate the workspace */
  AllocateWSpace(ctrl, 10*graph->nvtxs);


  /* Partition and Remap */
  STARTTIMER(ctrl, ctrl->TotalTmr);

  ctrl->ipc_factor = *ipc2redist;
  ctrl->CoarsenTo  = gk_min(graph->gnvtxs+1,
      (gk_max(npes, *nparts) > 256 ? 20 : 50)*(*ncon)*gk_max(npes, *nparts));

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


/*************************************************************************/
/*! This function is the driver for the adaptive refinement mode of 
    ParMETIS 
*/
/*************************************************************************/
void Adaptive_Partition(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i;
  idx_t tewgt, tvsize;
  real_t gtewgt, gtvsize;
  real_t ubavg, lbavg, *lbvec;

  WCOREPUSH;

  lbvec = rwspacemalloc(ctrl, graph->ncon);

  /************************************/
  /* Set up important data structures */
  /************************************/
  CommSetup(ctrl, graph);

  ubavg   = ravg(graph->ncon, ctrl->ubvec);
  tewgt   = isum(graph->nedges, graph->adjwgt, 1);
  tvsize  = isum(graph->nvtxs, graph->vsize, 1);
  gtewgt  = (real_t) GlobalSESum(ctrl, tewgt) + 1.0/graph->gnvtxs;  /* The +1/graph->gnvtxs were added to remove any FPE */
  gtvsize = (real_t) GlobalSESum(ctrl, tvsize) + 1.0/graph->gnvtxs;
  ctrl->redist_factor = ctrl->redist_base * ((gtewgt/gtvsize)/ ctrl->edge_size_ratio);

  IFSET(ctrl->dbglvl, DBG_PROGRESS, rprintf(ctrl, "[%6"PRIDX" %8"PRIDX" %5"PRIDX" %5"PRIDX"][%"PRIDX"]\n", 
        graph->gnvtxs, GlobalSESum(ctrl, graph->nedges), GlobalSEMin(ctrl, graph->nvtxs), GlobalSEMax(ctrl, graph->nvtxs), ctrl->CoarsenTo));

  if (graph->gnvtxs < 1.3*ctrl->CoarsenTo ||
     (graph->finer != NULL && graph->gnvtxs > graph->finer->gnvtxs*COARSEN_FRACTION)) {

    AllocateRefinementWorkSpace(ctrl, 2*graph->nedges);

    /***********************************************/
    /* Balance the partition on the coarsest graph */
    /***********************************************/
    graph->where = ismalloc(graph->nvtxs+graph->nrecv, -1, "graph->where");
    icopy(graph->nvtxs, graph->home, graph->where);

    ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
    lbavg = ravg(graph->ncon, lbvec);

    if (lbavg > ubavg + 0.035 && ctrl->partType != REFINE_PARTITION)
      Balance_Partition(ctrl, graph);

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

    /* check if no coarsening took place */
    if (graph->finer == NULL) {
      ComputePartitionParams(ctrl, graph);
      KWayBalance(ctrl, graph, graph->ncon);
      KWayAdaptiveRefine(ctrl, graph, NGR_PASSES);
    }
  }
  else {
    /*******************************/
    /* Coarsen it and partition it */
    /*******************************/
    switch (ctrl->ps_relation) {
      case PARMETIS_PSR_COUPLED:
        Match_Local(ctrl, graph);
        break;
      case PARMETIS_PSR_UNCOUPLED:
      default:
        Match_Global(ctrl, graph);
        break;
    }

    Adaptive_Partition(ctrl, graph->coarser);

    /********************************/
    /* project partition and refine */
    /********************************/
    ProjectPartition(ctrl, graph);
    ComputePartitionParams(ctrl, graph);

    if (graph->ncon > 1 && graph->level < 4) {
      ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
      lbavg = ravg(graph->ncon, lbvec);

      if (lbavg > ubavg + 0.025) {
        KWayBalance(ctrl, graph, graph->ncon);
      }
    }

    KWayAdaptiveRefine(ctrl, graph, NGR_PASSES);

    if (ctrl->dbglvl&DBG_PROGRESS) {
      ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
      rprintf(ctrl, "nvtxs: %10"PRIDX", cut: %8"PRIDX", balance: ", 
          graph->gnvtxs, graph->mincut);
      for (i=0; i<graph->ncon; i++) 
        rprintf(ctrl, "%.3"PRREAL" ", lbvec[i]);
      rprintf(ctrl, "\n");
    }
  }

  WCOREPOP;
}

