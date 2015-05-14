/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * stat.c
 *
 * This file computes various statistics
 *
 * Started 7/25/97
 * George
 *
 * $Id: stat.c 10578 2011-07-14 18:10:15Z karypis $
 *
 */

#include <parmetislib.h>



/*************************************************************************
* This function computes the balance of the partitioning
**************************************************************************/
void ComputeSerialBalance(ctrl_t *ctrl, graph_t *graph, idx_t *where, real_t *ubvec)
{
  idx_t i, j, nvtxs, ncon, nparts;
  idx_t *pwgts, *tvwgts, *vwgt;
  real_t *tpwgts, maximb;

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  vwgt   = graph->vwgt;
  nparts = ctrl->nparts;
  tpwgts = ctrl->tpwgts;

  pwgts = ismalloc(nparts*ncon, 0, "pwgts");
  tvwgts = ismalloc(ncon, 0, "tvwgts");

  for (i=0; i<graph->nvtxs; i++) {
    for (j=0; j<ncon; j++) {
      pwgts[where[i]*ncon+j] += vwgt[i*ncon+j];
      tvwgts[j] += vwgt[i*ncon+j];
    }
  }

  /* The +1 in the following code is to deal with bad cases of tpwgts[i*ncon+j] == 0 */
  for (j=0; j<ncon; j++) {
    maximb = 0.0;
    for (i=0; i<nparts; i++)
      maximb =gk_max(maximb, (1.0+(real_t)pwgts[i*ncon+j])/(1.0+(tpwgts[i*ncon+j]*(real_t)tvwgts[j])));
    ubvec[j] = maximb;
  }

  gk_free((void **)&pwgts, (void **)&tvwgts, LTERM);
}


/*************************************************************************
* This function computes the balance of the partitioning
**************************************************************************/
void ComputeParallelBalance(ctrl_t *ctrl, graph_t *graph, idx_t *where, real_t *ubvec)
{
  idx_t i, j, nvtxs, ncon, nparts;
  real_t *nvwgt, *lnpwgts, *gnpwgts, *lminvwgts, *gminvwgts;
  real_t *tpwgts, maximb;

  WCOREPUSH;

  ncon   = graph->ncon;
  nvtxs  = graph->nvtxs;
  nvwgt  = graph->nvwgt;
  nparts = ctrl->nparts;
  tpwgts = ctrl->tpwgts;

  lminvwgts = rset(ncon, 1.0, rwspacemalloc(ctrl, ncon)); 
  gminvwgts = rwspacemalloc(ctrl, ncon); 
  lnpwgts   = rset(nparts*ncon, 0.0, rwspacemalloc(ctrl, nparts*ncon));
  gnpwgts   = rwspacemalloc(ctrl, nparts*ncon);

  for (i=0; i<nvtxs; i++) {
    for (j=0; j<ncon; j++) {
      lnpwgts[where[i]*ncon+j] += nvwgt[i*ncon+j];

      /* The following is to deal with tpwgts[] that are 0.0 for certain partitions/constraints */
      lminvwgts[j] = (nvwgt[i*ncon+j] > 0.0 && lminvwgts[j] > nvwgt[i*ncon+j] ? nvwgt[i*ncon+j] : lminvwgts[j]);
    }
  }

  gkMPI_Allreduce((void *)(lnpwgts), (void *)(gnpwgts), nparts*ncon, REAL_T, MPI_SUM, ctrl->comm);
  gkMPI_Allreduce((void *)(lminvwgts), (void *)(gminvwgts), ncon, REAL_T, MPI_MIN, ctrl->comm);

  /* The +gminvwgts[j] in the following code is to deal with bad cases of tpwgts[i*ncon+j] == 0 */
  for (j=0; j<ncon; j++) {
    maximb = 0.0;
    for (i=0; i<nparts; i++)
      maximb =gk_max(maximb, (gminvwgts[j]+gnpwgts[i*ncon+j])/(gminvwgts[j]+tpwgts[i*ncon+j]));
    ubvec[j] = maximb;
  }

  WCOREPOP;
}


/*************************************************************************
* This function prints a matrix
**************************************************************************/
void Mc_PrintThrottleMatrix(ctrl_t *ctrl, graph_t *graph, real_t *matrix)
{
  idx_t i, j;

  for (i=0; i<ctrl->npes; i++) {
    if (i == ctrl->mype) {
      for (j=0; j<ctrl->npes; j++)
        printf("%.3"PRREAL" ", matrix[j]);
      printf("\n");
      fflush(stdout);
    }
    gkMPI_Barrier(ctrl->comm);
  }

  if (ctrl->mype == 0) {
    printf("****************************\n");
    fflush(stdout);
  }
  gkMPI_Barrier(ctrl->comm);

  return;
}


/***********************************************************************************/
/*! This function prints post-partitioning information 
   */
/***********************************************************************************/
void PrintPostPartInfo(ctrl_t *ctrl, graph_t *graph, idx_t movestats)
{
  idx_t i, j, ncon, nmoved, maxin, maxout, nparts;
  real_t maximb, *tpwgts;

  ncon = graph->ncon;

  nparts = ctrl->nparts;
  tpwgts = ctrl->tpwgts;

  rprintf(ctrl, "Final %3"PRIDX"-way Cut: %6"PRIDX" \tBalance: ", nparts, graph->mincut);

  for (j=0; j<ncon; j++) {
    for (maximb=0.0, i=0; i<nparts; i++) 
      maximb = gk_max(maximb, graph->gnpwgts[i*ncon+j]/tpwgts[i*ncon+j]);
    rprintf(ctrl, "%.3"PRREAL" ", maximb);
  }

  if (movestats) {
    Mc_ComputeMoveStatistics(ctrl, graph, &nmoved, &maxin, &maxout);
    rprintf(ctrl, "\nNMoved: %"PRIDX" %"PRIDX" %"PRIDX" %"PRIDX"\n", 
        nmoved, maxin, maxout, maxin+maxout);
  }
  else {
    rprintf(ctrl, "\n");
  }
}



/*************************************************************************
* This function computes movement statistics for adaptive refinement
* schemes
**************************************************************************/
void ComputeMoveStatistics(ctrl_t *ctrl, graph_t *graph, idx_t *nmoved, idx_t *maxin, idx_t *maxout)
{
  idx_t i, j, nvtxs;
  idx_t *vwgt, *where;
  idx_t *lpvtxs, *gpvtxs;

  nvtxs = graph->nvtxs;
  vwgt = graph->vwgt;
  where = graph->where;

  lpvtxs = ismalloc(ctrl->nparts, 0, "ComputeMoveStatistics: lpvtxs");
  gpvtxs = ismalloc(ctrl->nparts, 0, "ComputeMoveStatistics: gpvtxs");

  for (j=i=0; i<nvtxs; i++) {
    lpvtxs[where[i]]++;
    if (where[i] != ctrl->mype)
      j++;
  }

  /* PrintVector(ctrl, ctrl->npes, 0, lpvtxs, "Lpvtxs: "); */

  gkMPI_Allreduce((void *)lpvtxs, (void *)gpvtxs, ctrl->nparts, IDX_T, MPI_SUM, ctrl->comm);

  *nmoved = GlobalSESum(ctrl, j);
  *maxout = GlobalSEMax(ctrl, j);
  *maxin = GlobalSEMax(ctrl, gpvtxs[ctrl->mype]-(nvtxs-j));

  gk_free((void **)&lpvtxs, (void **)&gpvtxs, LTERM);
}


