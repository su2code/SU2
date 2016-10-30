/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * balancemylink.c
 *
 * This file contains code that implements the edge-based FM refinement
 *
 * Started 7/23/97
 * George
 *
 * $Id: balancemylink.c 10542 2011-07-11 16:56:22Z karypis $
 */

#include <parmetislib.h>
#define	PE	0

/*************************************************************************
* This function performs an edge-based FM refinement
**************************************************************************/
idx_t BalanceMyLink(ctrl_t *ctrl, graph_t *graph, idx_t *home, idx_t me,
          idx_t you, real_t *flows, real_t maxdiff, real_t *diff_cost, 
          real_t *diff_lbavg, real_t avgvwgt)
{
  idx_t h, i, ii, j, k, mype;
  idx_t nvtxs, ncon;
  idx_t nqueues, minval, maxval, higain, vtx, edge, totalv;
  idx_t from, to, qnum, index, nchanges, cut, tmp;
  idx_t pass, nswaps, nmoves, multiplier;
  idx_t *xadj, *vsize, *adjncy, *adjwgt, *where, *ed, *id;
  idx_t *hval, *nvpq, *inq, *map, *rmap, *ptr, *myqueue, *changes;
  real_t *nvwgt, *lbvec, *pwgts, *tpwgts, *my_wgt;
  real_t newgain;
  real_t lbavg, bestflow, mycost;
  real_t ipc_factor, redist_factor, ftmp;
  rpq_t **queues;

  WCOREPUSH;

  gkMPI_Comm_rank(MPI_COMM_WORLD, &mype);

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj;
  nvwgt  = graph->nvwgt;
  vsize  = graph->vsize;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  ipc_factor    = ctrl->ipc_factor;
  redist_factor = ctrl->redist_factor;

  hval    = iwspacemalloc(ctrl, nvtxs);
  id      = iwspacemalloc(ctrl, nvtxs);
  ed      = iwspacemalloc(ctrl, nvtxs);
  map     = iwspacemalloc(ctrl, nvtxs);
  rmap    = iwspacemalloc(ctrl, nvtxs);
  myqueue = iwspacemalloc(ctrl, nvtxs);
  changes = iwspacemalloc(ctrl, nvtxs);

  lbvec  = rwspacemalloc(ctrl, ncon);
  pwgts  = rset(2*ncon, 0.0, rwspacemalloc(ctrl, 2*ncon));
  tpwgts = rwspacemalloc(ctrl, 2*ncon);
  my_wgt = rset(ncon, 0.0, rwspacemalloc(ctrl, ncon));

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

  /* we don't want any tpwgts to be less than zero */
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

  /*******************************/
  /* insert vertices into queues */
  /*******************************/
  minval = maxval = 0;
  multiplier = 1;
  for (i=0; i<ncon; i++) {
    multiplier *= (i+1);
    maxval += i*multiplier;
    minval += (ncon-1-i)*multiplier;
  }

  nqueues = maxval-minval+1;
  nvpq    = iset(nqueues, 0, iwspacemalloc(ctrl, nqueues));
  ptr     = iwspacemalloc(ctrl, nqueues+1);
  inq     = iwspacemalloc(ctrl, 2*nqueues);
  queues  = (rpq_t **)(wspacemalloc(ctrl, sizeof(rpq_t *)*2*nqueues));

  for (i=0; i<nvtxs; i++)
    hval[i] = Mc_HashVwgts(ctrl, ncon, nvwgt+i*ncon) - minval;

  for (i=0; i<nvtxs; i++)
    nvpq[hval[i]]++;

  for (ptr[0]=0, i=0; i<nqueues; i++)
    ptr[i+1] = ptr[i] + nvpq[i];

  for (i=0; i<nvtxs; i++) {
    map[i] = ptr[hval[i]];
    rmap[ptr[hval[i]]++] = i;
  }

  SHIFTCSR(i, nqueues, ptr);


  /* initialize queues */
  for (i=0; i<nqueues; i++)
    if (nvpq[i] > 0) {
      queues[i]         = rpqCreate(nvpq[i]);
      queues[nqueues+i] = rpqCreate(nvpq[i]);
    }

  /* compute internal/external degrees */
  iset(nvtxs, 0, id);
  iset(nvtxs, 0, ed);
  for (j=0; j<nvtxs; j++) {
    for (k=xadj[j]; k<xadj[j+1]; k++) {
      if (where[adjncy[k]] == where[j])
        id[j] += adjwgt[k];
      else 
        ed[j] += adjwgt[k];
    }
  }

  nswaps = 0;
  for (pass=0; pass<N_MOC_BAL_PASSES; pass++) {
    iset(nvtxs, -1, myqueue); 
    iset(nqueues*2, 0, inq);

    /* insert vertices into correct queues */
    for (j=0; j<nvtxs; j++) {
      index = (where[j] == me) ? 0 : nqueues;

      newgain = ipc_factor*(real_t)(ed[j]-id[j]);
      if (home[j] == me || home[j] == you) {
        if (where[j] == home[j])
          newgain -= redist_factor*(real_t)vsize[j];
        else
          newgain += redist_factor*(real_t)vsize[j];
      }

      rpqInsert(queues[hval[j]+index], map[j]-ptr[hval[j]], newgain);
      myqueue[j] = (where[j] == me) ? 0 : 1;
      inq[hval[j]+index]++;
    }

    /* bestflow = rfavg(ncon, flows); */
    for (j=0, h=0; h<ncon; h++) {
      if (fabs(flows[h]) > fabs(flows[j])) 
        j = h;
    }
    bestflow = fabs(flows[j]);

    nchanges = nmoves = 0;
    for (ii=0; ii<nvtxs/2; ii++) {
      from = -1;
      Mc_DynamicSelectQueue(ctrl, nqueues, ncon, me, you, inq, flows, 
          &from, &qnum, minval, avgvwgt, maxdiff);

      /* can't find a vertex in one subdomain, try the other */
      if (from != -1 && qnum == -1) {
        from = (from == me) ? you : me;

        if (from == me) {
          for (j=0; j<ncon; j++) {
            if (flows[j] > avgvwgt)
              break;
          }
        }
        else {
          for (j=0; j<ncon; j++) {
            if (flows[j] < -1.0*avgvwgt)
              break;
          }
        }

        if (j != ncon)
          Mc_DynamicSelectQueue(ctrl, nqueues, ncon, me, you, inq, flows, 
              &from, &qnum, minval, avgvwgt, maxdiff);
      }

      if (qnum == -1)
        break;

      to = (from == me) ? you : me;
      index = (from == me) ? 0 : nqueues;
      higain = rpqGetTop(queues[qnum+index]);
      inq[qnum+index]--;
      ASSERT(higain != -1);

      /*****************/
      /* make the swap */
      /*****************/
      vtx = rmap[higain+ptr[qnum]];
      myqueue[vtx] = -1;
      where[vtx] = to;
      nswaps++;
      nmoves++;

      /* update the flows */
      for (j=0; j<ncon; j++)
        flows[j] += (to == me) ? nvwgt[vtx*ncon+j] : -1.0*nvwgt[vtx*ncon+j];
 
      /* ftmp = rfavg(ncon, flows); */
      for (j=0, h=0; h<ncon; h++) {
        if (fabs(flows[h]) > fabs(flows[j])) 
          j = h;
      }
      ftmp = fabs(flows[j]);

      if (ftmp < bestflow) {
        bestflow = ftmp;
        nchanges = 0;
      }
      else {
        changes[nchanges++] = vtx;
      }

      gk_SWAP(id[vtx], ed[vtx], tmp);

      for (j=xadj[vtx]; j<xadj[vtx+1]; j++) {
        edge = adjncy[j];

        tmp = (to == where[edge] ? adjwgt[j] : -adjwgt[j]);
        INC_DEC(id[edge], ed[edge], tmp);

        if (myqueue[edge] != -1) {
          newgain = ipc_factor*(real_t)(ed[edge]-id[edge]);
          if (home[edge] == me || home[edge] == you) {
            if (where[edge] == home[edge])
              newgain -= redist_factor*(real_t)vsize[edge];
            else
              newgain += redist_factor*(real_t)vsize[edge];
          }

          rpqUpdate(queues[hval[edge]+(nqueues*myqueue[edge])], 
              map[edge]-ptr[hval[edge]], newgain);
        }
      }
    }

    /****************************/
    /* now go back to best flow */
    /****************************/
    nswaps -= nchanges;
    nmoves -= nchanges;
    for (i=0; i<nchanges; i++) {
      vtx = changes[i];
      from = where[vtx];
      where[vtx] = to = (from == me) ? you : me;

     gk_SWAP(id[vtx], ed[vtx], tmp);
      for (j=xadj[vtx]; j<xadj[vtx+1]; j++) {
        edge = adjncy[j];
        tmp = (to == where[edge] ? adjwgt[j] : -adjwgt[j]);
        INC_DEC(id[edge], ed[edge], tmp);
      }
    }

    for (i=0; i<nqueues; i++) {
      if (nvpq[i] > 0) {
        rpqReset(queues[i]);
        rpqReset(queues[i+nqueues]);
      }
    }

    if (nmoves == 0)
      break;
  }

  /***************************/
  /* compute 2-way imbalance */
  /***************************/
  for (i=0; i<nvtxs; i++) {
    if (where[i] == me) {
      for (h=0; h<ncon; h++)
        my_wgt[h] += nvwgt[i*ncon+h];
    }
  }

  for (i=0; i<ncon; i++) {
    ftmp =  (pwgts[i]+pwgts[ncon+i])/2.0;
    if (ftmp != 0.0)
      lbvec[i] = fabs(my_wgt[i]-tpwgts[i]) / ftmp;
    else
      lbvec[i] = 0.0;
  }
  lbavg = ravg(ncon, lbvec);
  *diff_lbavg = lbavg;

  /****************/
  /* compute cost */
  /****************/
  cut = totalv = 0;
  for (i=0; i<nvtxs; i++) {
    if (where[i] != home[i])
      totalv += vsize[i];

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (where[adjncy[j]] != where[i])
          cut += adjwgt[j];
      }
  }
  cut /= 2;
  mycost = cut*ipc_factor + totalv*redist_factor;
  *diff_cost = mycost;

  /* free memory */
  for (i=0; i<nqueues; i++) {
    if (nvpq[i] > 0) {
      rpqDestroy(queues[i]);
      rpqDestroy(queues[i+nqueues]);
    }
  }

  WCOREPOP;

  return nswaps;
}

