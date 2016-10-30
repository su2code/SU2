/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * redomylink.c
 *
 * This file contains code that implements the edge-based FM refinement
 *
 * Started 7/23/97
 * George
 *
 * $Id: redomylink.c 10542 2011-07-11 16:56:22Z karypis $
 */

#include <parmetislib.h>
#define	PE	0

/*************************************************************************
* This function performs an edge-based FM refinement
**************************************************************************/
void RedoMyLink(ctrl_t *ctrl, graph_t *graph, idx_t *home, idx_t me,
         idx_t you, real_t *flows, real_t *sr_cost, real_t *sr_lbavg)
{
  idx_t h, i, r;
  idx_t nvtxs, nedges, ncon;
  idx_t  pass, lastseed, totalv;
  idx_t *xadj, *adjncy, *adjwgt, *where, *vsize;
  idx_t *costwhere, *lbwhere, *selectwhere;
  idx_t *ed, *id, *bndptr, *bndind, *perm;
  real_t *nvwgt, mycost;
  real_t lbavg, *lbvec;
  real_t best_lbavg, other_lbavg = -1.0, bestcost, othercost = -1.0;
  real_t *npwgts, *pwgts, *tpwgts;
  real_t ipc_factor, redist_factor, ftmp;
  idx_t mype;
  gkMPI_Comm_rank(MPI_COMM_WORLD, &mype);


  WCOREPUSH;

  nvtxs  = graph->nvtxs;
  nedges = graph->nedges;
  ncon   = graph->ncon;
  xadj   = graph->xadj;
  nvwgt  = graph->nvwgt;
  vsize  = graph->vsize;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  ipc_factor    = ctrl->ipc_factor;
  redist_factor = ctrl->redist_factor;

  /* set up data structures */
  id     = graph->sendind = iwspacemalloc(ctrl, nvtxs);
  ed     = graph->recvind = iwspacemalloc(ctrl, nvtxs);
  bndptr = graph->sendptr = iwspacemalloc(ctrl, nvtxs);
  bndind = graph->recvptr = iwspacemalloc(ctrl, nvtxs);

  costwhere = iwspacemalloc(ctrl, nvtxs);
  lbwhere   = iwspacemalloc(ctrl, nvtxs);
  perm      = iwspacemalloc(ctrl, nvtxs);

  lbvec  = rwspacemalloc(ctrl, ncon);
  pwgts  = rset(2*ncon, 0.0, rwspacemalloc(ctrl, 2*ncon));
  npwgts = rwspacemalloc(ctrl, 2*ncon);
  tpwgts = rwspacemalloc(ctrl, 2*ncon);

  graph->gnpwgts = npwgts;

  RandomPermute(nvtxs, perm, 1);
  icopy(nvtxs, where, costwhere);
  icopy(nvtxs, where, lbwhere);

  /* compute target pwgts */
  for (h=0; h<ncon; h++) {
    tpwgts[h]      = -1.0*flows[h];
    tpwgts[ncon+h] = flows[h];
  }

  for (i=0; i<nvtxs; i++) {
    if (where[i] == me) {
      for (h=0; h<ncon; h++) {
        tpwgts[h] += nvwgt[i*ncon+h];
        pwgts[h]  += nvwgt[i*ncon+h];
      }
    }
    else {
      ASSERT(where[i] == you);
      for (h=0; h<ncon; h++) {
        tpwgts[ncon+h] += nvwgt[i*ncon+h];
        pwgts[ncon+h]  += nvwgt[i*ncon+h];
      }
    }
  }

  /* we don't want any weights to be less than zero */
  for (h=0; h<ncon; h++) {
    if (tpwgts[h] < 0.0) {
      tpwgts[ncon+h] += tpwgts[h];
      tpwgts[h] = 0.0;
    }

    if (tpwgts[ncon+h] < 0.0) {
      tpwgts[h] += tpwgts[ncon+h];
      tpwgts[ncon+h] = 0.0;
    }
  } 


  /* now compute new bisection */
  bestcost = (real_t)isum(nedges, adjwgt, 1)*ipc_factor + 
             (real_t)isum(nvtxs, vsize, 1)*redist_factor;
  best_lbavg = 10.0;

  lastseed = 0;
  for (pass=N_MOC_REDO_PASSES; pass>0; pass--) {
    iset(nvtxs, 1, where);

    /* find seed vertices */
    r = perm[lastseed] % nvtxs;
    lastseed = (lastseed+1) % nvtxs;
    where[r] = 0;

    Mc_Serial_Compute2WayPartitionParams(ctrl, graph);
    Mc_Serial_Init2WayBalance(ctrl, graph, tpwgts);
    Mc_Serial_FM_2WayRefine(ctrl, graph, tpwgts, 4);
    Mc_Serial_Balance2Way(ctrl, graph, tpwgts, 1.02);
    Mc_Serial_FM_2WayRefine(ctrl, graph, tpwgts, 4);

    for (i=0; i<nvtxs; i++)
      where[i] = (where[i] == 0) ? me : you;

    for (i=0; i<ncon; i++) {
      ftmp = (pwgts[i]+pwgts[ncon+i])/2.0;
      if (ftmp != 0.0)
        lbvec[i] = fabs(npwgts[i]-tpwgts[i])/ftmp;
      else
        lbvec[i] = 0.0;
    }
    lbavg = ravg(ncon, lbvec);

    totalv = 0;
    for (i=0; i<nvtxs; i++)
      if (where[i] != home[i])
        totalv += vsize[i];

    mycost = (real_t)(graph->mincut)*ipc_factor + (real_t)totalv*redist_factor;

    if (bestcost >= mycost) {
      bestcost = mycost;
      other_lbavg = lbavg;
      icopy(nvtxs, where, costwhere);
    }

    if (best_lbavg >= lbavg) {
      best_lbavg = lbavg;
      othercost = mycost;
      icopy(nvtxs, where, lbwhere);
    }
  }

  if (other_lbavg <= .05) {
    selectwhere = costwhere;
    *sr_cost = bestcost;
    *sr_lbavg = other_lbavg;
  }
  else {
    selectwhere = lbwhere;
    *sr_cost = othercost;
    *sr_lbavg = best_lbavg;
  }

  icopy(nvtxs, selectwhere, where);

  WCOREPOP;
}

