/*!
 * Copyright 1997, Regents of the University of Minnesota
 *
 * \file
 * \brief This is the entry point of parallel ordering routines
 *
 * \date Started 8/1/2008
 * \author George Karypis
 * \version\verbatim $Id: ometis.c 10666 2011-08-04 05:22:36Z karypis $ \endverbatime
 *
 */

#include <parmetislib.h>


/***********************************************************************************/
/*! This function is the entry point of the parallel ordering algorithm. 
    It simply translates the arguments to the tunable version. */
/***********************************************************************************/
int ParMETIS_V3_NodeND(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, 
         idx_t *numflag, idx_t *options, idx_t *order, idx_t *sizes,
         MPI_Comm *comm)
{
  idx_t status;
  idx_t seed   = (options != NULL && options[0] != 0 ? options[PMV3_OPTION_SEED] : -1);
  idx_t dbglvl = (options != NULL && options[0] != 0 ? options[PMV3_OPTION_DBGLVL] : -1);

  /* Check the input parameters and return if an error */
  status = CheckInputsNodeND(vtxdist, xadj, adjncy, numflag, options, order, sizes, comm);
  if (GlobalSEMinComm(*comm, status) == 0) 
    return METIS_ERROR;


  ParMETIS_V32_NodeND(vtxdist, xadj, adjncy, 
      /*vwgt=*/NULL, 
      numflag, 
      /*mtype=*/NULL, 
      /*rtype=*/NULL, 
      /*p_nseps=*/NULL, 
      /*s_nseps=*/NULL, 
      /*ubfrac=*/NULL,
      /*seed=*/(options==NULL || options[0] == 0 ? NULL : &seed),
      /*dbglvl=*/(options==NULL || options[0] == 0 ? NULL : &dbglvl),
      order, sizes, comm);

  return METIS_OK;
}


/***********************************************************************************/
/*! This function is the entry point of the tunable parallel ordering algorithm. 
    This is the main ordering algorithm and implements a multilevel nested 
    dissection ordering approach. 
*/
/***********************************************************************************/
int ParMETIS_V32_NodeND(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
              idx_t *numflag, idx_t *mtype, idx_t *rtype, idx_t *p_nseps, idx_t *s_nseps, 
              real_t *ubfrac, idx_t *seed, idx_t *idbglvl, idx_t *order, idx_t *sizes, 
              MPI_Comm *comm)
{
  idx_t i, npes, mype, dbglvl, status, wgtflag=0;
  ctrl_t *ctrl;
  graph_t *graph, *mgraph;
  idx_t *morder;
  size_t curmem;

  gkMPI_Comm_size(*comm, &npes);
  gkMPI_Comm_rank(*comm, &mype);

  /* Deal with poor vertex distributions */
  if (GlobalSEMinComm(*comm, vtxdist[mype+1]-vtxdist[mype]) < 1) {
    printf("Error: Poor vertex distribution (processor with no vertices).\n");
    return METIS_ERROR;
  }

  status = METIS_OK;
  gk_malloc_init();
  curmem = gk_GetCurMemoryUsed();

  /* Setup the ctrl */
  ctrl = SetupCtrl(PARMETIS_OP_KMETIS, NULL, 1, 5*npes, NULL, NULL, *comm);

  dbglvl = (idbglvl == NULL ? 0 : *idbglvl);

  ctrl->dbglvl = dbglvl;
  STARTTIMER(ctrl, ctrl->TotalTmr);
  ctrl->dbglvl = 0;
  /*=======================================================================*/
  /*! Compute the initial k-way partitioning */
  /*=======================================================================*/
  /* Setup the graph */
  if (*numflag > 0) 
    ChangeNumbering(vtxdist, xadj, adjncy, order, npes, mype, 1);

  graph = SetupGraph(ctrl, 1, vtxdist, xadj, NULL, NULL, adjncy, NULL, 0);

  /* Allocate workspace */
  AllocateWSpace(ctrl, 10*graph->nvtxs);


  /* Compute the partitioning */
  ctrl->CoarsenTo = gk_min(vtxdist[npes]+1, 200*gk_max(npes, ctrl->nparts));
  if (seed != NULL) 
    ctrl->seed = (*seed == 0 ? mype : (*seed)*mype);

  Global_Partition(ctrl, graph);

  /* Collapse the number of partitions to be from 0..npes-1 */
  for (i=0; i<graph->nvtxs; i++)
    graph->where[i] = graph->where[i]%npes;
  ctrl->nparts = npes;

  /* Put back the real vertex weights */
  if (vwgt) {
    gk_free((void **)&graph->vwgt, LTERM);
    graph->vwgt      = vwgt;
    graph->free_vwgt = 0;
    wgtflag = 2;
  }


  /*=======================================================================*/
  /*! Move the graph according to the partitioning */
  /*=======================================================================*/
  STARTTIMER(ctrl, ctrl->MoveTmr);

  mgraph = MoveGraph(ctrl, graph);

  /* compute nvwgts for the moved graph */
  SetupGraph_nvwgts(ctrl, mgraph);

  STOPTIMER(ctrl, ctrl->MoveTmr);


  /*=======================================================================*/
  /*! Now compute an ordering of the moved graph */
  /*=======================================================================*/
  ctrl->optype    = PARMETIS_OP_OMETIS;
  ctrl->partType  = ORDER_PARTITION;
  ctrl->mtype     = (mtype  == NULL ? PARMETIS_MTYPE_GLOBAL  : *mtype);
  ctrl->rtype     = (rtype  == NULL ? PARMETIS_SRTYPE_2PHASE : *rtype);
  ctrl->p_nseps   = (p_nseps  == NULL ? 1 : *p_nseps);
  ctrl->s_nseps   = (s_nseps  == NULL ? 1 : *s_nseps);
  ctrl->ubfrac    = (ubfrac == NULL ? ORDER_UNBALANCE_FRACTION : *ubfrac);
  ctrl->dbglvl    = dbglvl;
  ctrl->ipart     = ISEP_NODE;
  ctrl->CoarsenTo = gk_min(graph->gnvtxs-1,
                       gk_max(1500*npes, graph->gnvtxs/(5*NUM_INIT_MSECTIONS*npes)));

  morder = imalloc(mgraph->nvtxs, "ParMETIS_NodeND: morder");
  MultilevelOrder(ctrl, mgraph, morder, sizes);

  /* Invert the ordering back to the original graph */
  ProjectInfoBack(ctrl, graph, order, morder);

  STOPTIMER(ctrl, ctrl->TotalTmr);
  IFSET(dbglvl, DBG_TIME, PrintTimingInfo(ctrl));
  IFSET(dbglvl, DBG_TIME, gkMPI_Barrier(ctrl->gcomm));

  gk_free((void **)&morder, LTERM);
  FreeGraph(mgraph);
  FreeInitialGraphAndRemap(graph);

  /* If required, restore the graph numbering */
  if (*numflag > 0) 
    ChangeNumbering(vtxdist, xadj, adjncy, order, npes, mype, 0);

  goto DONE;  /* Remove the warning for now */


DONE:
  FreeCtrl(&ctrl);
  if (gk_GetCurMemoryUsed() - curmem > 0) {
    printf("ParMETIS appears to have a memory leak of %zdbytes. Report this.\n",
        (ssize_t)(gk_GetCurMemoryUsed() - curmem));
  }
  gk_malloc_cleanup(0);

  return (int)status;
}


/*********************************************************************************/
/*!
  This is the top level ordering routine. 
  \param order is the computed ordering.
  \param sizes is the 2*nparts array that will store the sizes of each subdomains 
               and the sizes of the separators at each level. Note that the 
               top-level separator is stores at \c sizes[2*nparts-2].
*/
/*********************************************************************************/
void MultilevelOrder(ctrl_t *ctrl, graph_t *graph, idx_t *order, idx_t *sizes)
{
  idx_t i, nparts, nvtxs, npes;
  idx_t *perm, *lastnode, *morder, *porder;
  graph_t *mgraph;

  nvtxs = graph->nvtxs;

  npes = 1<<log2Int(ctrl->npes); /* # of nested dissection levels = floor(log_2(npes)) */

  perm     = imalloc(nvtxs, "MultilevelOrder: perm");
  lastnode = ismalloc(4*npes, -1, "MultilevelOrder: lastnode");

  for (i=0; i<nvtxs; i++) 
    perm[i] = i;
  lastnode[2] = graph->gnvtxs;

  iset(nvtxs, -1, order);

  /* This is used as a pointer to the end of the sizes[] array (i.e., >=nparts)
     that has not yet been filled in so that the separator sizes of the succesive
     levels will be stored correctly. It is used in LabelSeparatos() */
  sizes[0] = 2*npes-1;

  graph->where = ismalloc(nvtxs, 0, "MultilevelOrder: graph->where");

  for (nparts=2; nparts<=npes; nparts*=2) {
    ctrl->nparts = nparts;

    Order_Partition_Multiple(ctrl, graph);

    LabelSeparators(ctrl, graph, lastnode, perm, order, sizes);

    CompactGraph(ctrl, graph, perm);

    if (ctrl->CoarsenTo < 100*nparts) {
      ctrl->CoarsenTo = 1.5*ctrl->CoarsenTo;
    }
    ctrl->CoarsenTo = gk_min(ctrl->CoarsenTo, graph->gnvtxs-1);
  }


  /*-----------------------------------------------------------------
   / Move the graph so that each processor gets its partition 
   -----------------------------------------------------------------*/
  IFSET(ctrl->dbglvl, DBG_TIME, gkMPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MoveTmr));

  CommSetup(ctrl, graph);
  graph->ncon = 1; /* needed for MoveGraph */
  mgraph = MoveGraph(ctrl, graph);

  /* Fill in the sizes[] array for the local part. Just the vtxdist of the mgraph */
  for (i=0; i<npes; i++)
    sizes[i] = mgraph->vtxdist[i+1]-mgraph->vtxdist[i];

  porder = imalloc(graph->nvtxs, "MultilevelOrder: porder");
  morder = imalloc(mgraph->nvtxs, "MultilevelOrder: morder");

  IFSET(ctrl->dbglvl, DBG_TIME, gkMPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MoveTmr));

  /* Find the local ordering */
  if (ctrl->mype < npes)
    LocalNDOrder(ctrl, mgraph, morder, lastnode[2*(npes+ctrl->mype)]-mgraph->nvtxs);

  /* Project the ordering back to the before-move graph */
  ProjectInfoBack(ctrl, graph, porder, morder);

  /* Copy the ordering from porder to order using perm */
  for (i=0; i<graph->nvtxs; i++) {
    PASSERT(ctrl, order[perm[i]] == -1);
    order[perm[i]] = porder[i];
  }


  FreeGraph(mgraph);
  gk_free((void **)&perm, (void **)&lastnode, (void **)&porder, (void **)&morder, LTERM);

  /* PrintVector(ctrl, 2*npes-1, 0, sizes, "SIZES"); */
}


/**************************************************************************/
/*! This is the top-level driver of the multiple multisection ordering
    code. */
/***************************************************************************/
void Order_Partition_Multiple(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, sid, iter, nvtxs, nparts, nlevels;
  idx_t *xadj, *adjncy, *where, *gpwgts, *imap;
  idx_t *bestseps, *bestwhere, *origwhere;

  CommSetup(ctrl, graph);

  nparts = ctrl->nparts;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;

  bestseps  = ismalloc(2*nparts, -1, "Order_Partition_Multiple: bestseps");
  bestwhere = imalloc(nvtxs+graph->nrecv, "Order_Partition_Multiple: bestwhere");

  origwhere = graph->where;

  for (nlevels=-1, iter=0; iter<ctrl->p_nseps; iter++) {
    graph->where = imalloc(nvtxs, "Order_Partition_Multiple: where");
    icopy(nvtxs, origwhere, graph->where);

    Order_Partition(ctrl, graph, &nlevels, 0);

    where  = graph->where;
    gpwgts = graph->gpwgts;
    /* Update the where[] vectors of the subdomains that improved */
    for (i=0; i<nvtxs; i++) {
      sid = (where[i] < nparts ? nparts + where[i] - (where[i]%2) : where[i]);
      if (iter == 0 || bestseps[sid] > gpwgts[sid])
        bestwhere[i] = where[i];
    }
    /* Update the size of the separators for those improved subdomains */
    for (i=0; i<nparts; i+=2) {
      sid = nparts+i;
      if (iter == 0 || bestseps[sid] > gpwgts[sid]) 
        bestseps[sid] = gpwgts[sid];
    }

    /* free all the memory allocated for coarsening/refinement, but keep
       the setup fields so that they will not be re-computed */
    FreeNonGraphNonSetupFields(graph);
  }

  graph->where = bestwhere;
  AllocateNodePartitionParams(ctrl, graph);
  ComputeNodePartitionParams(ctrl, graph);

  for (i=0; i<nparts; i+=2) 
    PASSERT(ctrl, bestseps[nparts+i] == graph->gpwgts[nparts+i]);

  gk_free((void **)&bestseps, &origwhere, LTERM);

  /* PrintVector(ctrl, 2*nparts-1, 0, bestseps, "bestseps"); */

}


/**************************************************************************/
/*! The driver of the multilvelel separator finding algorithm */
/**************************************************************************/
void Order_Partition(ctrl_t *ctrl, graph_t *graph, idx_t *nlevels, idx_t clevel)
{

  CommSetup(ctrl, graph);
  graph->ncon = 1;

  IFSET(ctrl->dbglvl, DBG_PROGRESS, rprintf(ctrl, "[%6"PRIDX" %8"PRIDX" %5"PRIDX" %5"PRIDX"][%"PRIDX"][%"PRIDX"]\n",
        graph->gnvtxs, GlobalSESum(ctrl, graph->nedges), GlobalSEMin(ctrl, graph->nvtxs),
        GlobalSEMax(ctrl, graph->nvtxs), ctrl->CoarsenTo,
        GlobalSEMax(ctrl, imax(graph->nvtxs, graph->vwgt))));

  if ((*nlevels != -1 && *nlevels == clevel) ||
      (*nlevels == -1 && 
       ((graph->gnvtxs < 1.66*ctrl->CoarsenTo) || 
        (graph->finer != NULL && graph->gnvtxs > graph->finer->gnvtxs*COARSEN_FRACTION)))) {
    /* set the nlevels to where coarsening stopped */
    *nlevels = clevel;

    /* Compute the initial npart-way multisection */
    InitMultisection(ctrl, graph);

    if (graph->finer == NULL) { /* Do that only if no-coarsening took place */
      AllocateNodePartitionParams(ctrl, graph);
      ComputeNodePartitionParams(ctrl, graph);
      switch (ctrl->rtype) {
        case PARMETIS_SRTYPE_GREEDY:
          KWayNodeRefine_Greedy(ctrl, graph, NGR_PASSES, ctrl->ubfrac);
          break;
        case PARMETIS_SRTYPE_2PHASE:
          KWayNodeRefine2Phase(ctrl, graph, NGR_PASSES, ctrl->ubfrac);
          break;
        default:
          errexit("Unknown rtype of %"PRIDX"\n", ctrl->rtype);
      }
    }
  }
  else { /* Coarsen it and then partition it */
    switch (ctrl->mtype) {
      case PARMETIS_MTYPE_LOCAL:
        Match_Local(ctrl, graph);
        break;
      case PARMETIS_MTYPE_GLOBAL:
        Match_Global(ctrl, graph);
        break;
      default:
        errexit("Unknown mtype of %"PRIDX"\n", ctrl->mtype);
    }

    Order_Partition(ctrl, graph->coarser, nlevels, clevel+1);

    ProjectPartition(ctrl, graph);
    AllocateNodePartitionParams(ctrl, graph);
    ComputeNodePartitionParams(ctrl, graph);

    switch (ctrl->rtype) {
      case PARMETIS_SRTYPE_GREEDY:
        KWayNodeRefine_Greedy(ctrl, graph, NGR_PASSES, ctrl->ubfrac);
        break;
      case PARMETIS_SRTYPE_2PHASE:
        KWayNodeRefine2Phase(ctrl, graph, NGR_PASSES, ctrl->ubfrac);
        break;
      default:
        errexit("Unknown rtype of %"PRIDX"\n", ctrl->rtype);
    }
  }
}



/*********************************************************************************/
/*! This function is used to assign labels to the nodes in the separators. 
    It uses the appropriate entry in the lastnode array to select label boundaries
    and adjusts it for the next level. */
/*********************************************************************************/
void LabelSeparators(ctrl_t *ctrl, graph_t *graph, idx_t *lastnode, idx_t *perm, 
         idx_t *order, idx_t *sizes) 
{ 
  idx_t i, nvtxs, nparts, sid; idx_t *where, *lpwgts, *gpwgts, *sizescan;

  nparts = ctrl->nparts;

  nvtxs  = graph->nvtxs;
  where  = graph->where;
  lpwgts = graph->lpwgts;
  gpwgts = graph->gpwgts;

  if (ctrl->dbglvl&DBG_INFO) { 
    if (ctrl->mype == 0) {
      printf("SepWgts:  ");
      for (i=0; i<nparts; i+=2)
        printf(" %"PRIDX" [%"PRIDX" %"PRIDX"]", gpwgts[nparts+i], gpwgts[i], gpwgts[i+1]);
      printf("\n");
    }
    gkMPI_Barrier(ctrl->comm);
  }

  /* Compute the local size of the separator. This is required in case the 
     graph has vertex weights */
  iset(2*nparts, 0, lpwgts);
  for (i=0; i<nvtxs; i++) 
    lpwgts[where[i]]++;

  sizescan = imalloc(2*nparts, "LabelSeparators: sizescan");

  /* Perform a Prefix scan of the separator sizes to determine the boundaries */
  gkMPI_Scan((void *)lpwgts, (void *)sizescan, 2*nparts, IDX_T, MPI_SUM, ctrl->comm);
  gkMPI_Allreduce((void *)lpwgts, (void *)gpwgts, 2*nparts, IDX_T, MPI_SUM, ctrl->comm);

#ifdef DEBUG_ORDER
  PrintVector(ctrl, 2*nparts, 0, lpwgts, "Lpwgts");
  PrintVector(ctrl, 2*nparts, 0, sizescan, "SizeScan");
  PrintVector(ctrl, 2*nparts, 0, lastnode, "LastNode");
#endif

  /* Fillin the sizes[] array. See the comment on MultilevelOrder() on the 
     purpose of the sizes[0] value. */
  for (i=nparts-2; i>=0; i-=2) 
    sizes[--sizes[0]] = gpwgts[nparts+i];

  if (ctrl->dbglvl&DBG_INFO) { 
    if (ctrl->mype == 0) {
      printf("SepSizes: ");
      for (i=0; i<nparts; i+=2)
        printf(" %"PRIDX" [%"PRIDX" %"PRIDX"]", gpwgts[nparts+i], gpwgts[i], gpwgts[i+1]);
      printf("\n");
    }
    gkMPI_Barrier(ctrl->comm);
  }

  for (i=0; i<2*nparts; i++)
    sizescan[i] -= lpwgts[i];

  /* Assign the order[] values to the separator nodes */
  for (i=0; i<nvtxs; i++) {
    if (where[i] >= nparts) {
      sid = where[i];
      sizescan[sid]++;
      PASSERT(ctrl, order[perm[i]] == -1);
      order[perm[i]] = lastnode[sid] - sizescan[sid];
      /*myprintf(ctrl, "order[%"PRIDX"] = %"PRIDX", %"PRIDX"\n", perm[i], order[perm[i]], sid); */
    }
  }

  /* Update lastnode array */
  icopy(2*nparts, lastnode, sizescan);
  for (i=0; i<nparts; i+=2) {
    lastnode[2*nparts+2*i]     = sizescan[nparts+i]-gpwgts[nparts+i]-gpwgts[i+1];
    lastnode[2*nparts+2*(i+1)] = sizescan[nparts+i]-gpwgts[nparts+i];
    /*myprintf(ctrl, "lastnode: %"PRIDX" %"PRIDX"\n", lastnode[2*nparts+2*i], * lastnode[2*nparts+2*(i+1)]);*/
  }

  gk_free((void **)&sizescan, LTERM);

}



/*************************************************************************
* This function compacts a graph by removing the vertex separator
**************************************************************************/
void CompactGraph(ctrl_t *ctrl, graph_t *graph, idx_t *perm)
{
  idx_t i, j, l, nvtxs, cnvtxs, cfirstvtx, nparts, npes; 
  idx_t *xadj, *adjncy, *adjwgt, *vtxdist, *where;
  idx_t *cmap, *cvtxdist, *newwhere;

  WCOREPUSH;

  nparts = ctrl->nparts;
  npes   = ctrl->npes;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  if (graph->cmap == NULL)
    graph->cmap = imalloc(nvtxs+graph->nrecv, "CompactGraph: cmap");
  cmap = graph->cmap;

  vtxdist = graph->vtxdist;

  /*************************************************************
  * Construct the cvtxdist of the contracted graph. Uses the fact
  * that lpwgts stores the local non separator vertices.
  **************************************************************/
  cnvtxs   = isum(nparts, graph->lpwgts, 1);
  cvtxdist = iwspacemalloc(ctrl, npes+1);

  gkMPI_Allgather((void *)&cnvtxs, 1, IDX_T, (void *)cvtxdist, 1, IDX_T, 
      ctrl->comm);
  MAKECSR(i, npes, cvtxdist);

#ifdef DEBUG_ORDER
  PrintVector(ctrl, npes+1, 0, cvtxdist, "cvtxdist");
#endif


  /*************************************************************
  * Construct the cmap vector 
  **************************************************************/
  cfirstvtx = cvtxdist[ctrl->mype];

  /* Create the cmap of what you know so far locally */
  for (cnvtxs=0, i=0; i<nvtxs; i++) {
    if (where[i] < nparts) {
      perm[cnvtxs] = perm[i];
      cmap[i] = cfirstvtx + cnvtxs++;
    }
  }

  CommInterfaceData(ctrl, graph, cmap, cmap+nvtxs);


  /*************************************************************
  * Finally, compact the graph
  **************************************************************/
  newwhere = imalloc(cnvtxs, "CompactGraph: newwhere");
  cnvtxs = l = 0;
  for (i=0; i<nvtxs; i++) {
    if (where[i] < nparts) {
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        PASSERT(ctrl, where[i] == where[adjncy[j]] || where[adjncy[j]] >= nparts);
        if (where[i] == where[adjncy[j]]) {
          adjncy[l]   = cmap[adjncy[j]];
          adjwgt[l++] = adjwgt[j];
        }
      }

      xadj[cnvtxs] = l;
      graph->vwgt[cnvtxs] = graph->vwgt[i];
      newwhere[cnvtxs]    = where[i];
      cnvtxs++;
    }
  }
  SHIFTCSR(i, cnvtxs, xadj);

  gk_free((void **)&graph->match, (void **)&graph->cmap, (void **)&graph->lperm, 
         (void **)&graph->where, (void **)&graph->label, (void **)&graph->ckrinfo,
         (void **)&graph->nrinfo, (void **)&graph->lpwgts, (void **)&graph->gpwgts, 
         (void **)&graph->sepind, (void **)&graph->peind,
         (void **)&graph->sendptr, (void **)&graph->sendind, 
         (void **)&graph->recvptr, (void **)&graph->recvind, 
         (void **)&graph->imap, (void **)&graph->rlens, (void **)&graph->slens, 
         (void **)&graph->rcand, (void **)&graph->pexadj, 
         (void **)&graph->peadjncy, (void **)&graph->peadjloc, LTERM);
 
  graph->nvtxs  = cnvtxs;
  graph->nedges = l;
  graph->gnvtxs = cvtxdist[npes];
  graph->where  = newwhere;
  icopy(npes+1, cvtxdist, graph->vtxdist);

  WCOREPOP;
}


/*************************************************************************/
/*! This function orders the locally stored graph using MMD. The vertices 
    will be ordered from firstnode onwards. */
/*************************************************************************/
void LocalNDOrder(ctrl_t *ctrl, graph_t *graph, idx_t *order, idx_t firstnode)
{
  idx_t i, j, nvtxs, firstvtx, lastvtx;
  idx_t *xadj, *adjncy;
  idx_t *perm, *iperm;
  idx_t numflag=0, options[METIS_NOPTIONS];

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->SerialTmr));
  WCOREPUSH;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;

  firstvtx = graph->vtxdist[ctrl->mype];
  lastvtx  = graph->vtxdist[ctrl->mype+1];

  /* Relabel the vertices so that they are in local index space */
  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      PASSERT(ctrl, adjncy[j]>=firstvtx && adjncy[j]<lastvtx);
      adjncy[j] -= firstvtx;
    }
  }

  perm  = iwspacemalloc(ctrl, nvtxs+5);
  iperm = iwspacemalloc(ctrl, nvtxs+5);

  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NSEPS] = ctrl->s_nseps;

  METIS_NodeND(&nvtxs, xadj, adjncy, graph->vwgt, options, perm, iperm);

  for (i=0; i<nvtxs; i++) {
    PASSERT(ctrl, iperm[i]>=0 && iperm[i]<nvtxs);
    order[i] = firstnode+iperm[i];
  }

  WCOREPOP;
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->SerialTmr));
}



